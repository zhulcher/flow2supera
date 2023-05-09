import h5py as h5
import numpy as np
from LarpixParser import event_parser as EventParser
from LarpixParser.util import detector_configuration

class InputEvent:
    event_id = -1
    mc_packets_assn = None
    tracks  = None
    packets = None
    trajectories = None
    t0 = -1
    first_track_id = -1
    event_separator = ''

class InputReader:
    
    def __init__(self,parser_run_config, input_files=None):
        self._mc_packets_assn = None
        self._packets = None
        self._tracks = None
        self._trajectories = None
        self._vertices = None
        self._packet2event = None
        self._event_ids = None
        self._event_t0s = None
        self._if_spill = False
        self._run_config = parser_run_config
        
        if input_files:
            self.ReadFile(input_files)
    

    def __len__(self):
        if self._event_ids is None: return 0
        return len(self._event_ids)


    def __iter__(self):
        for entry in range(len(self)):
            yield self.GetEntry(entry)


    def _correct_t0s(self,event_t0s,num_event):
        # compute dt.
        dt=event_t0s[1:]-event_t0s[:-1]
        print(f'    Found {(dt==0).sum()} duplicate T0 values (removing)' )
        print(f'    Entries removed: {np.where(dt==0)[0]+1}')
        # generate a mask for dt>0
        mask=np.insert(np.where(dt>0)[0]+1,0,0)
        # apply mask
        corrected_t0s = event_t0s[mask]
        return corrected_t0s

    
    def ReadFile(self,input_files,verbose=False):
        mc_packets_assn = []
        packets = []
        tracks  = []
        trajectories = []
        vertices = []
        
        if type(input_files) == str:
            input_files = [input_files]
            
        for f in input_files:
            with h5.File(f,'r') as fin:
                mc_packets_assn.append(fin['mc_packets_assn'][:])
                packets.append(fin['packets'][:])
                tracks.append(fin['tracks'][:])
                trajectories.append(fin['trajectories'][:])
                vertices.append(fin['vertices'][:])
                if verbose: print('Read-in:',f)
                
        self._mc_packets_assn = np.concatenate(mc_packets_assn)
        self._packets = np.concatenate(packets)
        self._tracks  = np.concatenate(tracks )
        self._trajectories = np.concatenate(trajectories)
        self._vertices = np.concatenate(vertices)
        
        # create mapping
        self._packet2event = EventParser.packet_to_eventid(self._mc_packets_assn, self._tracks)
        
        packet_mask = self._packet2event != -1
        ctr_packet  = len(self._packets)
        ctr_invalid_packet = ctr_packet - packet_mask.sum()
        if verbose:
            print('    %d (%.2f%%) packets without an event ID assignment. They will be ignored.' % (ctr_invalid_packet,
                                                                                                     ctr_invalid_packet/ctr_packet)
                 )
        
        # create a list of unique Event IDs
        self._event_ids = np.unique(self._packet2event[packet_mask]).astype(np.int64)
        if verbose:
            missing_ids = [i for i in np.arange(np.min(self._event_ids),np.max(self._event_ids)+1,1) if not i in self._event_ids]
            print('    %d unique event IDs found.' % len(self._event_ids))
            print('    Potentially missing %d event IDs %s' % (len(missing_ids),str(missing_ids)))
        
        # create a list of corresponding T0s        
        self._event_t0s = EventParser.get_t0_event(self._vertices,self._run_config)

        # Assert strong assumptions here
        # the number of readout should be same as the number of valid Event IDs
        if not len(self._event_ids) == self._event_t0s.shape[0]:
            raise ValueError(f'Mismatch in the number of unique Event IDs {len(self._event_ids)} and event T0 counts {self._event_t0s.shape[0]}')

        # Now it's safe to assume all readout groups for every event shares the same T0
        self._event_t0s = self._event_t0s.flatten()
        

    def GetEvent(self,event_id):
        
        index_loc = (self._event_ids == event_id).nonzero()[0]
        
        if len(index_loc) < 1:
            print('Event ID',event_id,'not found in the data')
            print('Invalid read request (returning None)')
            return None
        
        return GetEntry(index_loc[0])
        
    def GetEntry(self,index):
        
        if index >= len(self._event_ids):
            print('Entry',index,'is above allowed entry index (<%d)' % len(self._event_ids))
            print('Invalid read request (returning None)')
            return None
        
        # Now return event info for the found index
        result = InputEvent()

        result.event_separator = self._run_config['event_separator']
        
        result.event_id = self._event_ids[index]
        result.t0 = self._event_t0s[index]

        mask = self._packet2event == result.event_id
        
        result.packets = self._packets[mask]
        result.mc_packets_assn = self._mc_packets_assn[mask]
        
        mask = self._tracks[self._run_config['event_separator']] == result.event_id
        result.tracks = self._tracks[mask]
        
        result.first_track_id = mask.nonzero()[0][0]
        
        mask = self._trajectories[self._run_config['event_separator']] == result.event_id
        result.trajectories = self._trajectories[mask]
        
        return result  
