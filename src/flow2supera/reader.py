import h5py 
import h5flow
import numpy as np

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

class FlowReader:
    
    def __init__(self, parser_run_config, input_files=None):
        self._input_files = input_files
        if not isinstance(input_files, str):
            raise TypeError('Input file must be a str type')
        self._event_ids = None
        self._event_t0s = None
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
        t0s_path = '/combined/t0/'
        events_data_path = 'charge/events/data/'
        event_hit_indices_path = 'charge/events/ref/charge/calib_final_hits/ref_region/'
        calib_final_hits_path = 'charge/calib_final_hits/'
        calib_prompt_hits_path = 'charge/calib_prompt_hits/'
        backtracked_hits_path = 'mc_truth/calib_final_hit_backtrack/data'
        packets_path = 'charge/packets'
        interactions_path = 'mc_truth/interactions/data'
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
            self._event_ids = events_data['id']
            self._event_t0s = flow_manager[events_path, t0s_path]
            self._event_hit_indices = flow_manager[event_hit_indices_path]
            self._hits = flow_manager[calib_final_hits_path+'data']
            self._backtracked_hits = flow_manager[backtracked_hits_path]
            self._is_sim = 'mc_truth' in fin.keys()
            if self._is_sim:
                #self._segments = flow_manager[events_path,
                #                              calib_final_hits_path,
                #                              calib_prompt_hits_path,
                #                              packets_path,
                #                              segments_path]
                self._segments = flow_manager[segments_path+'data']
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
        max_contributors = 100
        hit_threshold = 0.0001
        #backtracked_hits = self._backtracked_hits
        # TODO Calculate the length of this in advance and use reserve; appending is slow!
        truth_dict = {
            'segment_ids': [],
            'trajectory_ids': [],
        }

        segment_ids = []
        trajectory_ids = []

        for i_bt, backtracked_hit in enumerate(backtracked_hits):
            for contrib in range(max_contributors):
                if abs(backtracked_hit['fraction'][contrib]) < hit_threshold: continue
                segment_id = backtracked_hit['segment_id'][contrib]
                segment = segments[segment_id]
                segment_ids.append(segment_id)
                trajectory_id = segment['traj_id']
                trajectory = trajectories[trajectory_id]
                trajectory_parent_id = trajectory['parent_id']
                trajectory_ids.append(trajectory_id)
                # Some trajectories' parents don't appear in the main trajectories
                # list, but need to be seen by the driver. Add them here explicitly.
                trajectory_ids.append(trajectory_parent_id)

        truth_dict['segment_ids'] = segment_ids
        # Trajectory IDs should increment in order
        truth_dict['trajectory_ids'] = sorted(trajectory_ids)

        return truth_dict

    def GetEvent(self, event_index):
        
        if event_index >= len(self._event_ids):
            print('Entry {} is above allowed entry index ({})'.format(event_index, len(self._event_ids)))
            print('Invalid read request (returning None)')
            return None
        
        result = InputEvent()

        result.event_id = self._event_ids[event_index]

        # t0s dtypes: ('id', 'ts', 'ts_err', 'type')
        # Use 'ts' for event timestamp
        result.t0 = self._event_t0s[result.event_id]['ts']

        result.hit_indices = self._event_hit_indices[result.event_id]
        hit_start_index = self._event_hit_indices[result.event_id][0]
        hit_stop_index  = self._event_hit_indices[result.event_id][1]
        result.hits = self._hits[hit_start_index:hit_stop_index]
        result.backtracked_hits = self._backtracked_hits[hit_start_index:hit_stop_index]

        truth_ids_dict = self.GetEventTruthFromHits(result.backtracked_hits, 
                                                    self._segments, 
                                                    self._trajectories)
        #result.trajectories = self._trajectories
        event_trajectory_ids = truth_ids_dict['trajectory_ids']
        trajectories_array = np.array(self._trajectories)
        result.trajectories = trajectories_array[np.isin(trajectories_array['traj_id'], event_trajectory_ids)]

        #result.segments = self._segments[result.event_id]
        event_segment_ids = truth_ids_dict['segment_ids']
        segments_array = np.array(self._segments)
        result.segments = segments_array[np.isin(segments_array['segment_id'], event_segment_ids)]

        #result.segments = self._segments
        #result.segments = self._segments[self._segments['event_id']==result.event_id]
        # Keep trajectories as-is and use the segments' traj_id to get event trajectories in driver
        #result.trajectories = self._trajectories[self._trajectories['event_id']==result.event_id]
        #result.trajectories = self._trajectories
        #result.interactions = self._interactions[result.event_id]
        #result.segment_index_min = mask.nonzero()[0][0]
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
        #print('trajectories:', input_event.trajectories)
        #print('segments:', input_event.segments)
        #print('trajectories:', input_event.trajectories)
        print('interactions in this event:', len(input_event.interactions))

