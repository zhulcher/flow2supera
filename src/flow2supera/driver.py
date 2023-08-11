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
        #max_trackid = max(data.trajectories['traj_id'].max(), data.segments['segment_id'].max())
        #self._trackid2idx.resize(int(max_trackid+1), supera.kINVALID_INDEX)
        for traj in data.trajectories:
            part_input = supera.ParticleInput()

            part_input.valid = True
            part_input.part  = self.TrajectoryToParticle(traj)
            part_input.part.id = supera_event.size()
            if self.GetLogger().verbose():
                if verbose:
                    print('  TrackID',part_input.part.trackid,
                          'PDG',part_input.part.pdg,
                          'Energy',part_input.part.energy_init)
            if traj['traj_id'] < 0:
                print('Negative track ID found',traj['traj_id'])
                raise ValueError
            #self._trackid2idx[int(traj['track_id'])] = part_input.part.id
            supera_event.push_back(part_input)
            
        if verbose:
            print("--- trajectory filling %s seconds ---" % (time.time() - start_time)) 
        start_time = time.time()  

        # 2. Fill parent information for ParticleInputs created in previous loop
        for i, part in enumerate(supera_event):
            traj = data.trajectories[i]

            parent=None            
            #if(part.part.parent_trackid < self._trackid2idx.size()):
            if(part.part.parent_trackid >= len(data.trajectories)): continue
            #parent_index = self._trackid2idx[part.part.parent_trackid]
            parent_index = event_trajectories[part.part.parent_trackid]
            if parent_index == supera.kINVALID_INDEX: continue

            parent = supera_event[parent_index].part
            part.part.parent_pdg = parent.pdg
                    
            self.SetProcessType(traj, part.part, parent)

        # TODO I think this should now loop over hits and backtracked hits using ref_region
        # 3. Loop over "voxels" (aka packets), get EDep from xyz and charge information,
        #    and store in pcloud
        #x, y, z, dE = HitParser.hit_parser_energy(data.t0, data.packets, self._geom_dict, self._run_config, switch_xz=True)

        start_time = time.time()

        # Loop over packets and decide particle trajectory segments that are associated with it.
        # Also calculate how much fraction of the packet value should be associated to this particle.

        # Flow files store up to 100 contributors per hit, but most of them are empty,
        # hence the threshold
        max_contributors = 100
        hit_threshold = 0.0001
        #hits = data.hits
        backtracked_hits = data.backtracked_hits
        for i_bt, backtracked_hit in enumerate(backtracked_hits):
            # TODO Loop over only non-zero contributors
            for contrib in range(max_contributors):
                if abs(backtracked_hit['fraction'][contrib]) < hit_threshold: continue

                reco_hit = data.hits[i_bt][contrib]

                # Record this packet
                #raw_edep.x, raw_edep.y, raw_edep.z = x[ip]*self._mm2cm, y[ip]*self._mm2cm, z[ip]*self._mm2cm
                #raw_edep.e = dE[ip]
                #self._edeps_all.push_back(raw_edep)
                #check_raw_sum += dE[ip]
                edep = supera.EDep()
                print('reco hit type:', type(reco_hit))
                print('reco hit shape:', reco_hit.shape)
                print('reco hit dtypes:', reco_hit.dtype.names)
                print('reco hit x:', reco_hit['x'])
                print('reco hit:', reco_hit)
                edep.x = reco_hit['x']
                edep.y = reco_hit['y']
                edep.z = reco_hit['z']
                edep.t = reco_hit['t_drift']

                #seg_idx = packet_segments[it]
                segment_id = backtracked_hit['segment_id'][contrib]
                print('segment id:', segment_id)
                print('data.segments len:', len(data.segments))
                segment = data.segments[segment_id]
                print('segment shape', segment.shape)
                print('segment dyptes', segment.dtype.names)
                print('segment dEdx', segment['dEdx'])
                print('segment', segment)
                edep.dedx = segment['dEdx']
                edep.e = reco_hit['E'] * backtracked_hit['fraction'][contrib]
                #supera_event[self._trackid2idx[int(seg['trackID'])]].pcloud.push_back(packet_edeps[it])
                supera_event[i_bt].pcloud.push_back(edep)

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
        p.id             = int(trajectory['event_id'])
        p.interaction_id = int(trajectory['vertex_id'])
        p.trackid        = int(trajectory['traj_id']) # TODO Is this rigth? What exactly is traj_id?
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
