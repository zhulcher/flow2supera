import h5py as h5
import h5flow
import numpy as np
#from LarpixParser import event_parser as EventParser
#from LarpixParser.util import detector_configuration

class InputEvent:
    event_id = -1
    #mc_packets_assn = None
    segments = None
    calib_final_hits  = None
    trajectories = None
    t0 = -1
    segment_index_min = -1
    event_separator = ''

class FlowReader:
    
    def __init__(self,parser_run_config, input_files=None):
        #self._mc_packets_assn = None
        self._calib_final_hits = None
        self._segments = None
        self._trajectories = None
        self._vertices = None
        self._packet2event = None
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
        #packets  = []
        calib_final_hits  = []
        segments = []
        trajectories = []
        t0s = []

        # H5Flow's H5FlowDataManager class associated datasets through references
        # These paths help us get the correct associations
        events_path = 'charge/events/'
        t0s_path = '/combined/t0/'
        calib_final_hits_path = 'charge/calib_final_hits/'
        calib_prompt_hits_path = 'charge/calib_prompt_hits/'
        packets_path = 'charge/packets'
        # TODO "tracks" will likely be renamed to "segments" soon in flow
        segments_path = '/mc_truth/tracks/'
        trajectories_path = '/mc_truth/trajectories/'

        if type(input_files) == str:
            input_files = [input_files]
        
        self._is_sim = True # TODO Should default to False
        for f in input_files:
            with h5.File(f, 'r') as fin:
                flow_manager = h5flow.data.H5FlowDataManager(fin, 'r')
                events = flow_manager[events_path+'data']
                calib_final_hits.append(flow_manager[events_path, 
                                                     calib_final_hits_path])
                t0s.append(flow_manager[events_path, t0s_path])
                #self._is_sim = 'mc_truth' in fin.keys()
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
                    
            # TODO Event IDs are not unique between files, I think. This needs to be 
            # checked. Maybe we just enforce one file per run of flow2supera? Seems
            # inefficient...
            self._events = np.concatenate(events)
            self._calib_final_hits = np.concatenate(calib_final_hits)

            if not self._is_sim:
                print('Currently only simulation is supoprted')
                raise NotImplementedError
            self._segments  = np.concatenate(segments )
            self._trajectories = np.concatenate(trajectories)
            self._event_t0s = np.concatenate(t0s)
            #self._vertices = np.concatenate(vertices)
        
    ### TODO ###
    # Fill event t0s
    # Fill event IDs

    def GetEvent(self,event_id):
        
        index_loc = (self._event_ids == event_id).nonzero()[0]
        
        if len(index_loc) < 1:
            print('Event ID',event_id,'not found in the data')
            print('Invalid read request (returning None)')
            return None
        
        return GetEntry(index_loc[0])

    def CheckIntegrity(self, data, ignore_bad_association=False):

        flag = True
        tid_range0 = np.array([t['trackID'] for t in data.trajectories])
        tid_range1 = np.array([s['trackID'] for s in data.segments    ])

        if tid_range0.max() < tid_range1.max():
            print('[ERROR] Max Track ID in the segments exceeds the maximum of the trajectories')
            flag = False

        if tid_range0.min() > tid_range1.min():
            print('[ERROR] Min Track ID in the segments is below the minimum of the trajectories')
            flag = False

        if not flag:
            return flag

        #seg_index = data.mc_packets_assn['track_ids']
        # check if max index is within the number of segments
        max_index = seg_index.max()
        min_index = seg_index[seg_index>-1].min()

        prefix = '[WARNING]' if ignore_bad_association else '[ERROR]'
        if min_index < data.segment_index_min:
            # Bad segment index on low end

            print(prefix,'Minimum segment index from the association:',min_index)
            print('        Index range of segments for this event:',data.segment_index_min,
                '=>',data.segment_index_min+len(data.segments))
            flag = False
            if ignore_bad_association:
                print('[WARNING] ignoring the bad association')
                seg_index[seg_index<data.segment_index_min] = -1
                #data.mc_packets_assn['track_ids'] = seg_index
                flag = True

        if (max_index - data.segment_index_min) >= len(data.segments):
            # Bad segment index on high end
            print(prefix,'Maximum segment index from the association:',max_index)
            print('        Index range of segments for this event:',data.segment_index_min,
                '=>',data.segment_index_min+len(data.segments))
            flag = False
            if ignore_bad_association:
                print('[WARNING] ignoring the bad association')
                seg_index[seg_index>=(data.segment_index_min+len(data.segments))] = -1
                #data.mc_packets_assn['track_ids'] = seg_index
                flag = True

        return flag


    def GetEntry(self,index):
        
        if index >= len(self._event_ids):
            print('Entry',index,'is above allowed entry index (<%d)' % len(self._event_ids))
            print('Invalid read request (returning None)')
            return None
        
        # Now return event info for the found index
        result = InputEvent()

        result.event_separator = self._run_config['event_separator']
        
        result.event_id = self._event_ids[index]
        result.t0 = self._event_t0s[result.event_id]

        mask = self._packet2event == result.event_id
        
        result.calib_final_hits = self._calib_final_hits[mask]
        #result.mc_packets_assn = self._mc_packets_assn[mask]
        
        mask = self._segments[self._run_config['event_separator']] == result.event_id
        result.segments = self._segments[mask]
        
        result.segment_index_min = mask.nonzero()[0][0]
        
        mask = self._trajectories[self._run_config['event_separator']] == result.event_id
        result.trajectories = self._trajectories[mask]
        
        return result  
