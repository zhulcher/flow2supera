import h5py 
import h5flow
import numpy as np

class InputEvent:
    event_id = -1
    segments = None
    hits_ref_region = None
    hits = None
    calib_final_hits  = None
    trajectories = None
    t0 = -1
    segment_index_min = -1
    event_separator = ''

class FlowReader:
    
    def __init__(self, parser_run_config, input_files=None):
        self._input_files = input_files
        if not isinstance(input_files, str):
            raise TypeError('Input file must be a str type')
        self._hits_ref_region = None
        self._hits = None
        self._backtracked_hits = None
        self._segments = None
        self._trajectories = None
        #self._vertices = None
        self._event_ids = None
        self._event_t0s = None
        self._interactions = None
        #self._if_spill = False
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
        hit_ref_region = []
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
        hit_ref_region_path = 'charge/events/ref/charge/calib_final_hits/ref_region/'
        t0s_path = '/combined/t0/'
        calib_final_hits_path = 'charge/calib_final_hits/'
        calib_prompt_hits_path = 'charge/calib_prompt_hits/'
        packets_path = 'charge/packets'
        interactions_path = 'mc_truth/interactions/'
        segments_path = 'mc_truth/segments/'
        trajectories_path = 'mc_truth/trajectories/data'

        #if type(input_files) == str:
        #    input_files = [input_files]
        
        self._is_sim = False 
        # TODO Currently only reading one input file at a time. Is it 
        # necessary to read multiple? If so, how to handle non-unique
        # event IDs?
        #for f in input_files:
        flow_manager = h5flow.data.H5FlowDataManager(input_files, 'r')
        with h5py.File(input_files, 'r') as fin:
            events = flow_manager[events_path]
            events_data = events['data']
            #event_ids.append(events_data['id'])
            self._event_ids = events_data['id']
            self._event_t0s = flow_manager[events_path, t0s_path]
            #self._hits_ref_region = fin[hit_ref_region_path][:]
            self._hits_ref_region = flow_manager[hit_ref_region_path]
            #calib_final_hits.append(flow_manager[events_path, 
            #                                     calib_final_hits_path])
            #self._hits = fin[calib_final_hits_path][:]
            self._hits = flow_manager[events_path, calib_final_hits_path]
            #t0s.append(flow_manager[events_path, t0s_path])
            self._is_sim = 'mc_truth' in fin.keys()
            if self._is_sim:
                #mc_packets_assn.append(fin['mc_packets_assn'][:])
                # TODO Truth backtracking should be streamlined starting in MiniRun4
                #segments.append(flow_manager[events_path,
                self._segments = flow_manager[events_path,
                                              calib_final_hits_path,
                                              calib_prompt_hits_path,
                                              packets_path,
                                              segments_path]
                #trajectories.append(flow_manager[events_path,
                #self._trajectories = flow_manager[events_path,
                #                                  calib_final_hits_path,
                #                                  calib_prompt_hits_path,
                #                                  packets_path,
                #                                  segments_path,
                #                                  trajectories_path]
                # Trajectories need to be separated manually since we want all
                # truth information, not just those mapped to reco
                self._trajectories = flow_manager[trajectories_path]

                self._interactions = flow_manager[interactions_path]
                
        # Stack datasets so that there's a "file index" preceding the event index
        #self._event_ids = np.stack(event_ids)
        #self._event_ids = np.concatenate(event_ids)
        #self._event_t0s = np.stack(t0s)
        #self._calib_final_hits = np.stack(calib_final_hits)
        #self._t0s = np.stack(t0s)
        #self._segments = np.stack(segments)
        #self._trajectories = np.stack(trajectories)

        print('len event_ids:', len(self._event_ids))
        print('len Reader:', len(self))

        if not self._is_sim:
            print('Currently only simulation is supoprted')
            raise NotImplementedError

    #def GetEntry(self, file_index, event_index):
    def GetEvent(self, event_index):
        
        if event_index >= len(self._event_ids):
            print('Entry {} is above allowed entry index ({})'.format(event_index, len(self._event_ids)))
            print('Invalid read request (returning None)')
            return None
        
        # Now return event info for the found index
        result = InputEvent()

        # TODO Is this necessary anymore?
        #result.event_separator = self._run_config['event_separator']
        
        result.event_id = self._event_ids[event_index]
        result.t0       = self._event_t0s[result.event_id]

        result.hits_ref_region = self._hits_ref_region[result.event_id]
        result.hits = self._hits[result.event_id]
        result.segments = self._segments[result.event_id]
        result.trajectories = self._trajectories[result.event_id]
        #result.mc_packets_assn = self._mc_packets_assn[mask]
        #result.segment_index_min = mask.nonzero()[0][0]
        
        return result  
