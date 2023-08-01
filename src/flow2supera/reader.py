import h5py 
import h5flow
import numpy as np

class InputEvent:
    event_id = -1
    segments = None
    calib_final_hits  = None
    trajectories = None
    t0 = -1
    segment_index_min = -1
    event_separator = ''

class FlowReader:
    
    def __init__(self, parser_run_config, input_files=None):
        self._calib_final_hits = None
        self._segments = None
        self._trajectories = None
        self._vertices = None
        self._event_ids = None
        self._event_t0s = None
        self._if_spill = False
        self._run_config = parser_run_config
        self._is_sim = False

        if input_files:
            self.ReadFile(input_files)

    def __len__(self):
        if self._event_ids is None: return 0
        return len(self._event_ids)

    def __iter__(self):
        for entry in range(len(self)):
            yield self.GetEntry(entry)

    #def _correct_t0s(self,event_t0s,num_event):
    #    # compute dt.
    #    dt=event_t0s[1:]-event_t0s[:-1]
    #    print(f'    Found {(dt==0).sum()} duplicate T0 values (removing)' )
    #    print(f'    Entries removed: {np.where(dt==0)[0]+1}')
    #    # generate a mask for dt>0
    #    mask=np.insert(np.where(dt>0)[0]+1,0,0)
    #    # apply mask
    #    corrected_t0s = event_t0s[mask]
    #    return corrected_t0s
    
    def ReadFile(self, input_files, verbose=False):
        event_ids = []
        calib_final_hits  = []
        segments = []
        trajectories = []
        t0s = []

        # H5Flow's H5FlowDataManager class associated datasets through references
        # These paths help us get the correct associations
        events_path = 'charge/events/'
        events_data_path = 'charge/events/data'
        t0s_path = '/combined/t0/'
        calib_final_hits_path = 'charge/calib_final_hits/'
        calib_prompt_hits_path = 'charge/calib_prompt_hits/'
        packets_path = 'charge/packets'
        # TODO "tracks" will likely be renamed to "segments" soon in flow
        segments_path = '/mc_truth/tracks/'
        trajectories_path = '/mc_truth/trajectories/'

        if type(input_files) == str:
            input_files = [input_files]
        
        self._is_sim = False 
        for f in input_files:
            flow_manager = h5flow.data.H5FlowDataManager(f, 'r')
            with h5py.File(f, 'r') as fin:
                events = flow_manager[events_path]
                events_data = events['data']
                event_ids.append(events_data['id'])
                calib_final_hits.append(flow_manager[events_path, 
                                                     calib_final_hits_path])
                t0s.append(flow_manager[events_path, t0s_path])
                self._is_sim = 'mc_truth' in fin.keys()
                if self._is_sim:
                    #mc_packets_assn.append(fin['mc_packets_assn'][:])
                    # TODO Truth backtracking should be streamlined starting in MiniRun4
                    segments.append(flow_manager[events_path,
                                                 calib_final_hits_path,
                                                 calib_prompt_hits_path,
                                                 packets_path,
                                                 segments_path])
                    trajectories.append(flow_manager[events_path,
                                                     calib_final_hits_path,
                                                     calib_prompt_hits_path,
                                                     packets_path,
                                                     segments_path,
                                                     trajectories_path])
                    if verbose: print('Read file:', fin )
                    
            # Stack datasets so that there's a "file index" preceding the event index
            self._event_ids = np.stack(event_ids)
            self._event_t0s = np.stack(t0s)
            self._calib_final_hits = np.stack(calib_final_hits)
            self._t0s = np.stack(t0s)
            self._segments = np.stack(segments)
            self._trajectories = np.stack(trajectories)

            if not self._is_sim:
                print('Currently only simulation is supoprted')
                raise NotImplementedError
        
    # TODO I think this isn't necessary anymore? 
    #def GetEvent(self, event_id):
    #    
    #    index_loc = (self._event_ids == event_id).nonzero()[0]
    #    
    #    if len(index_loc) < 1:
    #        print('Event ID',event_id,'not found in the data')
    #        print('Invalid read request (returning None)')
    #        return None
    #    
    #    return GetEntry(index_loc[0])

    def GetEntry(self, index):
        
        #event_id = -1
        #segments = None
        #calib_final_hits  = None
        #trajectories = None
        #t0 = -1
        #segment_index_min = -1
        #event_separator = ''

        if index >= len(self._event_ids):
            print('Entry',index,'is above allowed entry index (<%d)' % len(self._event_ids))
            print('Invalid read request (returning None)')
            return None
        
        # Now return event info for the found index
        result = InputEvent()

        # TODO Is this necessary anymore?
        result.event_separator = self._run_config['event_separator']
        
        result.event_id = self._event_ids[index]
        result.t0       = self._event_t0s[result.event_id]

        result.calib_final_hits = self._calib_final_hits[result.event_id]
        result.segments = self._segments[result.event_id]
        result.trajectories = self._trajectories[mask.event_id]
        #result.mc_packets_assn = self._mc_packets_assn[mask]
        #result.segment_index_min = mask.nonzero()[0][0]
        
        return result  
