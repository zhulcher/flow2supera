import h5py 
import h5flow
import numpy as np
from ROOT import supera
import cppyy

class InputEvent:
    event_id = -1
    segments = None
    hit_indices = None
    hits = None
    backtracked_hits = None
    calib_final_hits  = None
    trajectories = None
    interactions = None
    t0 = -1
    segment_index_min = -1
    event_separator = ''


class FlowChargeReader:
    
    def __init__(self, parser_run_config, input_files=None):
        self._input_files = input_files
        if not isinstance(input_files, str):
            raise TypeError('Input file must be a str type')
        self._event_ids = None
        self._event_t0s = None
        self._flash_t0s = None
        self._flash_ids = None
        self._event_hit_indices = None
        self._hits = None
        self._backtracked_hits = None
        self._segments = None
        self._trajectories = None
        self._interactions = None
        self._run_config = parser_run_config
        self._is_sim = False

        if input_files:
            self.ReadFile(input_files)

    def __len__(self):
        if self._event_ids is None: return 0
        return len(self._event_ids)

    
    def __iter__(self):
        for entry in range(len(self)):
            yield self.GetEvent(entry)

    def ReadFile(self, input_files, verbose=False):
        event_ids = []
        calib_final_hits  = []
        event_hit_indices = []
        hits = []
        backtracked_hits = []
        segments = []
        trajectories = []
        event_trajectories = []
        t0s = []

        print('Reading input file...')

        # H5Flow's H5FlowDataManager class associated datasets through references
        # These paths help us get the correct associations
        events_path = 'charge/events/'
        events_data_path = 'charge/events/data/'
        event_hit_indices_path = 'charge/events/ref/charge/calib_prompt_hits/ref_region/'
        calib_final_hits_path = 'charge/calib_final_hits/data'
        calib_prompt_hits_path = 'charge/calib_prompt_hits/data'
        backtracked_hits_path = 'mc_truth/calib_prompt_hit_backtrack/data'
        packets_path = 'charge/packets'
        interactions_path = 'mc_truth/interactions/data'
        segments_path = 'mc_truth/segments/data'
        trajectories_path = 'mc_truth/trajectories/data'

        self._is_sim = False 
        # TODO Currently only reading one input file at a time. Is it 
        # necessary to read multiple? If so, how to handle non-unique
        # event IDs?
        #for f in input_files:
        flow_manager = h5flow.data.H5FlowDataManager(input_files, 'r')
        with h5py.File(input_files, 'r') as fin:
            events = flow_manager[events_path]
            events_data = events['data']
            self._event_ids = events_data['id']
            self._event_t0s = events_data['ts_start']
            self._event_hit_indices = flow_manager[event_hit_indices_path]
            self._hits = flow_manager[calib_prompt_hits_path]
            self._backtracked_hits = flow_manager[backtracked_hits_path]
            self._is_sim = 'mc_truth' in fin.keys()
            if self._is_sim:
                #self._segments = flow_manager[events_path,
                #                              calib_final_hits_path,
                #                              calib_prompt_hits_path,
                #                              packets_path,
                #                              segments_path]
                self._segments = flow_manager[segments_path]
                self._trajectories = flow_manager[trajectories_path]
                self._interactions = flow_manager[interactions_path]

        # This next bit is only necessary if reading multiple files
        # Stack datasets so that there's a "file index" preceding the event index
        #self._event_ids = np.stack(event_ids)
        #self._event_ids = np.concatenate(event_ids)
        #self._event_t0s = np.stack(t0s)
        #self._calib_final_hits = np.stack(calib_final_hits)
        #self._t0s = np.stack(t0s)
        #self._segments = np.stack(segments)
        #self._trajectories = np.stack(trajectories)

        if not self._is_sim:
            print('Currently only simulation is supoprted')
            raise NotImplementedError

        
    # To truth associations go as hits -> segments -> trajectories
  
    def GetEventTruthFromHits(self, backtracked_hits, segments, trajectories):
        '''
        The Driver class needs to know the number of event trajectories in advance.
        This function uses the backtracked hits dataset to map hits->segments->trajectories
        and fills segment and trajectory IDs corresponding to hits. 
        '''
        truth_dict = {
            'segment_ids': [],
            'trajectory_ids': [],
        }

        trajectory_dict = {(traj['traj_id'], traj['event_id'], traj['vertex_id']): traj for traj in trajectories} #use all three to map the correct parent trajectory as 'parent_id' is the 'traj_id' which is not unique

        segment_ids = []
        trajectory_ids = []
        for backtracked_hit in backtracked_hits:
            for contrib in range(len(backtracked_hit['fraction'])):
                if abs(backtracked_hit['fraction'][contrib]) == 0: break

                segment_id = backtracked_hit['segment_id'][contrib]
                segment = segments[segment_id]
                segment_ids.append(segment_id)

                traj_id = segment['traj_id']
                event_id = segment['event_id']  
                vertex_id = segment['vertex_id']  

                trajectory_key = (traj_id, event_id, vertex_id)
                trajectory = trajectory_dict.get(trajectory_key)
                while trajectory is not None:
                    trajectory_ids.append(trajectory['file_traj_id'])
                    trajectory_parent_id = trajectory['parent_id'] 
                    if(trajectory_parent_id < 0): break #if <0, it is the parent
                    trajectory_parent_key = (trajectory_parent_id, event_id, vertex_id)
                    trajectory_parent = trajectory_dict.get(trajectory_parent_key)
                    trajectory = trajectory_parent
                # Some trajectories' parents don't appear in the main trajectories
                # list, but need to be seen by the driver. Add them here explicitly.
                
        truth_dict['segment_ids'] = segment_ids
        truth_dict['trajectory_ids'] = sorted(trajectory_ids)

        return truth_dict

    def GetEvent(self, event_index):
        
        if event_index >= len(self._event_ids):
            print('Entry {} is above allowed entry index ({})'.format(event_index, len(self._event_ids)))
            print('Invalid read request (returning None)')
            return None
        
        result = InputEvent()

        result.event_id = self._event_ids[event_index]


        result.t0 = self._event_t0s[result.event_id]*1e-7 #ts_start is in ticks and 0.1 microseconds per tick for charge readout

        result.hit_indices = self._event_hit_indices[result.event_id]
        hit_start_index = self._event_hit_indices[result.event_id][0]
        hit_stop_index  = self._event_hit_indices[result.event_id][1]
        result.hits = self._hits[hit_start_index:hit_stop_index]
        result.backtracked_hits = self._backtracked_hits[hit_start_index:hit_stop_index]

        truth_ids_dict = self.GetEventTruthFromHits(result.backtracked_hits, 
                                                    self._segments, 
                                                    self._trajectories)
        event_trajectory_ids = truth_ids_dict['trajectory_ids']
        trajectories_array = np.array(self._trajectories)
        result.trajectories = trajectories_array[np.isin(trajectories_array['file_traj_id'], event_trajectory_ids)]

        event_segment_ids = truth_ids_dict['segment_ids']
        segments_array = np.array(self._segments)
        result.segments = segments_array[np.isin(segments_array['segment_id'], event_segment_ids)]


        result.interactions = self._interactions
        
        return result  
 


    def EventDump(self, input_event):
        print('-----------EVENT DUMP-----------------')
        print('Event ID {}'.format(input_event.event_id))
        print('Event t0 {}'.format(input_event.t0))
        print('Event hit indices (start, stop):', input_event.hit_indices)
        print('Backtracked hits len:', len(input_event.backtracked_hits))
        print('hits shape:', input_event.hits.shape)
        print('segments in this event:', len(input_event.segments))
        print('trajectories in this event:', len(input_event.trajectories))
        print('interactions in this event:', len(input_event.interactions))

class FlowNeutrinoReader:
    
    def __init__(self, parser_run_config, input_files=None):
        self._input_files = input_files
        if not isinstance(input_files, str):
            raise TypeError('Input file must be a str type')

        self._ixns = None

       

        if input_files:
            self.ReadNeutrino(input_files)

    def __len__(self):
        if self._ixns is None: return 0
        return len(self._ixns)

    
    def __iter__(self):
        for entry in range(len(self)):
            yield self.GetNeutrino(entry)

    def ReadNeutrino(self, input_files, verbose=False):
        neutrino_info_path = 'mc_truth/interactions/data'
        flow_manager = h5flow.data.H5FlowDataManager(input_files, 'r')
        with h5py.File(input_files, 'r') as fin:
            self._ixns = flow_manager[neutrino_info_path]
    
    def GetNeutrino(self, event_index):
        
        if event_index >= len(self._ixns):
            print('Entry {} is above allowed entry index ({})'.format(event_index, len(self._ixns)))
            print('Invalid read request (returning None)')
            return None
        
        nu_result = supera.Neutrino()
        nu_result.id = event_index
        
        ixn = self._ixns[event_index]     
        nu_result.interaction_id = int(ixn['vertex_id']) 
        nu_result.target = int(ixn['target'])
        nu_result.vtx = supera.Vertex(ixn['vertex'][0], ixn['vertex'][1], ixn['vertex'][2], ixn['vertex'][3])
        nu_result.pdg_code = int(ixn['nu_pdg'])
        nu_result.lepton_pdg_code = int(ixn['lep_pdg'])  
        nu_result.energy_init = ixn['Enu']
        nu_result.theta = ixn['lep_ang']
        nu_result.momentum_transfer =  ixn['Q2']
        nu_result.momentum_transfer_mag =  ixn['q3']
        nu_result.energy_transfer =  ixn['q0']
        nu_result.bjorken_x = ixn['x']
        nu_result.inelasticity = ixn['y']
        nu_result.px = ixn['nu_4mom'][0]
        nu_result.py = ixn['nu_4mom'][1]       
        nu_result.pz = ixn['nu_4mom'][2]
        nu_result.lepton_p = ixn['lep_mom']
        if(ixn['isCC']): nu_result.current_type = 0
        else: nu_result.current_type = 1
        nu_result.interaction_mode = ixn['reaction']
            
        return nu_result  

