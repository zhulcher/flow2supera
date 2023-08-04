import time
import edep2supera, ROOT
from ROOT import supera, std, TG4TrajectoryPoint
import numpy as np
import LarpixParser
from LarpixParser import hit_parser as HitParser
import yaml
from yaml import Loader
import flow2supera

class SuperaDriver(edep2supera.edep2supera.SuperaDriver):

    LOG_KEYS = ('ass_saturation',
        'residual_q',
        'packet_total',
        'packet_ctr',
        'packet_noass',
        'packet_badass',
        'fraction_nan',
        'ass_frac')

    def __init__(self):
        super().__init__()
        self._geom_dict  = None
        self._run_config = None
        self._trackid2idx = std.vector('supera::Index_t')()
        self._allowed_detectors = std.vector('std::string')()
        self._edeps_unassociated = std.vector('supera::EDep')()
        self._edeps_all = std.vector('supera::EDep')()
        self._ass_distance_limit=0.4434*4.5
        self._ass_charge_limit=0.05
        self._log=None
        self._electron_energy_threshold=0
        self._estimate_pt_time=True
        self._ignore_bad_association=True
        print("Initialized SuperaDriver class")

    def parser_run_config(self):
        return self._run_config

    def log(self, data_holder):

        for key in self.LOG_KEYS:
            if key in data_holder:
                raise KeyError(f'Key {key} exists in the log data holder already.')
            data_holder[key]=[]
        self._log = data_holder

    def LoadPropertyConfigs(self, cfg_dict):

        # Expect only PropertyKeyword or (TileLayout,DetectorProperties). Not both.
        if cfg_dict.get('PropertyKeyword',None):
            print('hi')
            if cfg_dict.get('TileLayout',None) or cfg_dict.get('DetectorProperties',None):
                print('hi again')
                print('PropertyKeyword provided:', cfg_dict['PropertyKeyword'])
                print('But also founnd below:')
                for keyword in ['TileLayout','DetectorProperties']:
                    print('%s: "%s"' % (keyword,cfg_dict.get(keyword, None)))
                    print('Bool', bool(cfg_dict.get(keyword,None)))

                print('You cannot specify duplicated property infomration!')
                return False
            else:
                print('else')
                try:
                    print('try')
                    self._run_config, self._geom_dict = LarpixParser.util.detector_configuration(cfg_dict['PropertyKeyword'])
                except ValueError:
                    print('Failed to load with PropertyKeyword', cfg_dict['PropertyKeyword'])
                    print('Supported types:', LarpixParser.util.configuration_keywords())
                    return False
        else:
            print('PropertyKeyword missing. Loading TileLayout and DetectorProperties...')
            if not 'TileLayout' in cfg_dict:
                print('TileLayout not in the configuration data!')
                return False

            if not 'DetectorProperties' in cfg_dict:
                print('DetectorProperties not in the configuration data!')
                raise False

            self._geom_dict  = LarpixParser.util.load_geom_dict(cfg_dict['TileLayout'])
            self._run_config = LarpixParser.util.get_run_config(cfg_dict['DetectorProperties'])

        print('self._geom_dict', self._geom_dict)
        print('self._run_config', self._run_config)
        # Event separator default value needs to be set.
        # We repurpose "run_config" of EventParser to hold this attribute.
        self._run_config['event_separator'] = 'eventID'
        # Apply run config modification if requested
        run_config_mod = cfg_dict.get('ParserRunConfig',None)
        if run_config_mod:
            for key,val in run_config_mod.items():
                self._run_config[key]=val
        print('Returning property configs')
        return True

    def ConfigureFromFile(self,fname):
        with open(fname,'r') as f:
            cfg=yaml.load(f.read(),Loader=Loader)
            #if not self.LoadPropertyConfigs(cfg):
            #    raise ValueError('Failed to configure flow2supera!')
            self._electron_energy_threshold = cfg.get('ElectronEnergyThreshold',
                self._electron_energy_threshold
                )
            self._estimate_pt_time = cfg.get('EstimatePointTime',
                self._estimate_pt_time)
            self._ass_distance_limit = cfg.get('AssDistanceLimit',
                self._ass_distance_limit)
            self._ass_charge_limit = cfg.get('AssChargeLimit',
                self._ass_charge_limit)

        super().ConfigureFromFile(fname)

    def ConfigureFromText(self,txt):
        cfg=yaml.load(txt,Loader=Loader)
        if not self.LoadPropertyConfigs(cfg):
            raise ValueError('Failed to configure flow2supera!')
            self._electron_energy_threshold = cfg.get('ElectronEnergyThreshold',
                self._electron_energy_threshold
                )
        super().ConfigureFromText(txt)

    def ReadEvent(self, data, verbose=False):
        
        start_time = time.time()

        # initialize the new event record
        if not self._log is None:
            for key in self.LOG_KEYS:
                self._log[key].append(0)

        supera_event = supera.EventInput()
        supera_event.reserve(len(data.trajectories))
        
        # 1. Loop over trajectories, create one supera::ParticleInput for each
        #    store particle inputs in list to fill parent information later
        max_trackid = max(data.trajectories['trackID'].max(), data.segments['trackID'].max())
        self._trackid2idx.resize(int(max_trackid+1), supera.kINVALID_INDEX)
        for traj in data.trajectories:
            # print("traj",traj)
            part_input = supera.ParticleInput()

            part_input.valid = True
            part_input.part  = self.TrajectoryToParticle(traj)
            part_input.part.id = supera_event.size()
            #part_input.type  = self.GetG4CreationProcess(traj)
            #supera_event.append(part_input)
            if self.GetLogger().verbose():
                if verbose:
                    print('  TrackID',part_input.part.trackid,
                          'PDG',part_input.part.pdg,
                          'Energy',part_input.part.energy_init)
            if traj['trackID'] < 0:
                print('Negative track ID found',traj['trackID'])
                raise ValueError
            self._trackid2idx[int(traj['trackID'])] = part_input.part.id
            supera_event.push_back(part_input)
            
        if verbose:
            print("--- trajectory filling %s seconds ---" % (time.time() - start_time)) 
        start_time = time.time()  

        # 2. Fill parent information for ParticleInputs created in previous loop
        for i, part in enumerate(supera_event):
            traj = data.trajectories[i]

            parent=None            
            if(part.part.parent_trackid < self._trackid2idx.size()):
                parent_index = self._trackid2idx[part.part.parent_trackid]
                if not parent_index == supera.kINVALID_INDEX:
                    parent = supera_event[parent_index].part
                    part.part.parent_pdg = parent.pdg
                    
            self.SetProcessType(traj, part.part, parent)

        # TODO I think this is now loop over hits and backtracked hits using ref_region
        # 3. Loop over "voxels" (aka packets), get EDep from xyz and charge information,
        #    and store in pcloud
        #x, y, z, dE = HitParser.hit_parser_energy(data.t0, data.packets, self._geom_dict, self._run_config, switch_xz=True)

        start_time = time.time()

        ass_segments = data.mc_packets_assn['track_ids']
        ass_fractions = data.mc_packets_assn['fraction']
        
        # Loop over packets and decide particle trajectory segments that are associated with it.
        # Also calculate how much fraction of the packet value should be associated to this particle.
        check_raw_sum=0
        check_ana_sum=0
        for ip, packet in enumerate(data.packets):
            if verbose:
                print('*****************Packet', ip, '**********************')

            # If packet_type !=0 continue
            if packet['packet_type'] != 0: continue

            # Record this packet
            raw_edep = supera.EDep()
            raw_edep.x, raw_edep.y, raw_edep.z = x[ip]*self._mm2cm, y[ip]*self._mm2cm, z[ip]*self._mm2cm
            raw_edep.e = dE[ip]
            self._edeps_all.push_back(raw_edep)
            check_raw_sum += dE[ip]

            # Log the valid packet count
            if not self._log is None:
                self._log['packet_ctr'][-1]+=1

            # We analyze and modify segments and fractions, so make a copy
            packet_segments  = np.array(ass_segments[ip])
            packet_fractions = np.array(ass_fractions[ip])
            packet_edeps = [None] * len(packet_segments)

            # Initialize seg_flag once 
            if seg_flag is None:
                seg_flag = np.zeros(shape=(packet_segments.shape[0]),dtype=bool)
                seg_dist = np.zeros(shape=(packet_segments.shape[0]),dtype=float)
            seg_flag[:] = True
            seg_dist[:] = 1.e9

            if verbose:
                print('[INFO] Assessing packet',ip)
                print('       Segments :', packet_segments)
                print('       TrackIDs :', [data.segments[packet_segments[idx]]['trackID'] for idx in range(packet_segments.shape[0])])
                print('       Fractions:', ['%.3f' % f for f in packet_fractions])
                print('       Energy   : %.3f' % dE[ip])
                print('       Position :', ['%.3f' % f for f in [x[ip]*self._mm2cm,y[ip]*self._mm2cm,z[ip]*self._mm2cm]])
                print('       Distance :', ['%.3f' % f for f in seg_dist])

            # Post-Step1, check if there is any associated segments left.
            # Post-Step1.A: If none left, then register to the set of unassociated segments.
            if (seg_dist < self._ass_distance_limit).sum()<1:
                if verbose:
                    print('[WARNING] found a packet unassociated with any tracks')
                    print('          Packet %d Charge %.3f' % (ip,dE[ip]))
                edep = supera.EDep()
                edep.x,edep.y,edep.z,edep.e = x[ip]*self._mm2cm, y[ip]*self._mm2cm, z[ip]*self._mm2cm, dE[ip]
                self._edeps_unassociated.push_back(edep)
                if not self._log is None:
                    self._log['packet_noass'][-1] += 1
                check_ana_sum += edep.e
                continue

            # After this point, >=1 association will be found. Log if needed
            if not self._log is None:
                self._log['ass_frac'][-1] += 1

            # Post-Step1.B: If only one segment left, register to that packet.
            if (seg_dist < self._ass_distance_limit).sum()==1:
                if verbose:
                    print('[INFO] Registering the only segment within distance limit',packet_segments[it])
                it = np.argmin(seg_dist)
                seg_idx = packet_segments[it]
                seg = data.segments[seg_idx]
                packet_edeps[it].dedx = seg['dEdx']
                packet_edeps[it].e = dE[ip]
                supera_event[self._trackid2idx[int(seg['trackID'])]].pcloud.push_back(packet_edeps[it])
                check_ana_sum += packet_edeps[it].e
                continue

        return supera_event

    def TrajectoryToParticle(self, trajectory):
        ### What we have access to in new flow format: ###
        # ('event_id', 'vertex_id', 'traj_id', 'local_traj_id', 'parent_id', 'E_start', 'pxyz_start', 
        # 'xyz_start', 't_start', 'E_end', 'pxyz_end', 'xyz_end', 't_end', 'pdg_id', 
        # 'start_process', 'start_subprocess', 'end_process', 'end_subprocess')
        ###############################################################################################

        p = supera.Particle()
        # Larnd-sim stores a lot of these fields as numpy.uint32, 
        # but Supera/LArCV want a regular int, hence the type casting
        # TODO Is there a cleaner way to handle this?
        p.id             = int(trajectory['event_id'])
        p.interaction_id = trajectory['vertex_id']
        p.trackid        = int(trajectory['trackID'])
        p.pdg            = int(trajectory['pdg_id'])
        p.px = trajectory['pxyz_start'][0] 
        p.py = trajectory['pxyz_start'][1] 
        p.pz = trajectory['pxyz_start'][2]
        # TODO Verify this new thing works
        p.energy_init = trajectory['E_start']
        #p.energy_init = np.sqrt(pow(flow2supera.pdg2mass.pdg2mass(p.pdg),2) + 
        #                        pow(p.px,2) + pow(p.py,2) + pow(p.pz,2))
        # TODO Is this correct? Shouldn't the vertex be the interaction vertex?
        # And this should be p.start_pt or something?
        p.vtx    = supera.Vertex(trajectory['xyz_start'][0], 
                                 trajectory['xyz_start'][1], 
                                 trajectory['xyz_start'][2], 
                                 trajectory['t_start']
        )
        p.end_pt = supera.Vertex(trajectory['xyz_end'][0], 
                                 trajectory['xyz_end'][1],
                                 trajectory['xyz_end'][2], 
                                 trajectory['t_end']
        )

        traj_parent_id = trajectory['parent_id']
        # This now causes errors?
        #if traj_parent_id == -1: p.parent_trackid = supera.kINVALID_TRACKID
        if traj_parent_id == -1: p.parent_trackid = p.trackid
        else:                    p.parent_trackid = int(trajectory['parent_id'])

        if supera.kINVALID_TRACKID in [p.trackid, p.parent_trackid]:
            print('Unexpected to have an invalid track ID', p.trackid,
                  'or parent track ID', p.parent_trackid)
            raise ValueError
        
        return p
        
    def SetProcessType(self, edepsim_part, supera_part, supera_parent):
        pdg_code    = supera_part.pdg
        g4type_main = edepsim_part['start_process']
        g4type_sub  = edepsim_part['start_subprocess']
        
        supera_part.process = '%d::%d' % (int(g4type_main),int(g4type_sub))

        ke = np.sqrt(pow(supera_part.px,2)+pow(supera_part.py,2)+pow(supera_part.pz,2))

        dr = 1.e20
        if supera_parent:
            dx = (supera_parent.end_pt.pos.x - supera_part.vtx.pos.x)
            dy = (supera_parent.end_pt.pos.y - supera_part.vtx.pos.y)
            dz = (supera_parent.end_pt.pos.z - supera_part.vtx.pos.z)
            dr = dx + dy + dz
        
        if pdg_code == 2112 or pdg_code > 1000000000:
            supera_part.type = supera.kNeutron
        
        elif supera_part.trackid == supera_part.parent_trackid:
            supera_part.type = supera.kPrimary
            
        elif pdg_code == 22:
            supera_part.type = supera.kPhoton
        
        elif abs(pdg_code) == 11:
            
            if g4type_main == TG4TrajectoryPoint.G4ProcessType.kProcessElectromagetic:
                
                if g4type_sub == TG4TrajectoryPoint.G4ProcessSubtype.kSubtypeEMPhotoelectric:
                    supera_part.type = supera.kPhotoElectron
                    
                elif g4type_sub == TG4TrajectoryPoint.G4ProcessSubtype.kSubtypeEMComptonScattering:
                    supera_part.type = supera.kCompton
                
                elif g4type_sub == TG4TrajectoryPoint.G4ProcessSubtype.kSubtypeEMGammaConversion or \
                     g4type_sub == TG4TrajectoryPoint.G4ProcessSubtype.kSubtypeEMPairProdByCharged:
                    supera_part.type = supera.kConversion
                    
                elif g4type_sub == TG4TrajectoryPoint.G4ProcessSubtype.kSubtypeEMIonization:
                    
                    if abs(supera_part.parent_pdg) == 11:
                        supera_part.type = supera.kIonization
                        
                    elif abs(supera_part.parent_pdg) in [211,13,2212,321]:
                        supera_part.type = supera.kDelta

                    elif supera_part.parent_pdg == 22:
                        supera_part.type = supera.kCompton

                    else:
                        print("    WARNING: UNEXPECTED CASE for IONIZATION ")
                        print("      PDG",pdg_code,
                              "TrackId",edepsim_part['trackID'],
                              "Kinetic Energy",ke,
                              "Parent PDG",supera_part.parent_pdg ,
                              "Parent TrackId",edepsim_part['parentID'],
                              "G4ProcessType",g4type_main ,
                              "SubProcessType",g4type_sub)
                        supera_part.type = supera.kIonization
                #elif g4type_sub == 151:

                else:
                    print("    WARNING: UNEXPECTED EM SubType ")
                    print("      PDG",pdg_code,
                          "TrackId",edepsim_part['trackID'],
                          "Kinetic Energy",ke,
                          "Parent PDG",supera_part.parent_pdg ,
                          "Parent TrackId",edepsim_part['parentID'],
                          "G4ProcessType",g4type_main ,
                          "SubProcessType",g4type_sub)
                    raise ValueError
                    
            elif g4type_main == TG4TrajectoryPoint.G4ProcessType.kProcessDecay:
                #print("    WARNING: DECAY ")
                #print("      PDG",pdg_code,
                #      "TrackId",edepsim_part['trackID'],
                #      "Kinetic Energy",ke,
                #      "Parent PDG",supera_part.parent_pdg ,
                #      "Parent TrackId",edepsim_part['parentID'],
                #      "G4ProcessType",g4type_main ,
                #      "SubProcessType",g4type_sub)
                supera_part.type = supera.kDecay

            elif g4type_main == TG4TrajectoryPoint.G4ProcessType.kProcessHadronic and g4type_sub == 151 and dr<0.0001:
                if ke < self._electron_energy_threshold:
                    supera_part.type = supera.kIonization
                else:
                    supera_part.type = supera.kDecay
            
            else:
                print("    WARNING: Guessing the shower type as", "Compton" if ke < self._electron_energy_threshold else "OtherShower")
                print("      PDG",pdg_code,
                      "TrackId",edepsim_part['trackID'],
                      "Kinetic Energy",ke,
                      "Parent PDG",supera_part.parent_pdg ,
                      "Parent TrackId",edepsim_part['parentID'],
                      "G4ProcessType",g4type_main ,
                      "SubProcessType",g4type_sub)

                if ke < self._electron_energy_threshold:
                    supera_part.type = supera.kCompton
                else:
                    supera_part.type = supera.kOtherShower
        else:
            supera_part.type = supera.kTrack
