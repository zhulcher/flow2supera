import os
import numpy as np
import time
import flow2supera
import ROOT
from larcv import larcv
import yaml
from yaml import Loader
from edep2supera.utils import get_iomanager, larcv_meta, larcv_particle, larcv_neutrino

#from LarpixParser import event_parser as EventParser
from larcv import larcv

def get_iomanager(outname):
    import tempfile
    cfg='''                                                                                                                                          
IOManager: {                                                                                                                                         
  Verbosity:   2                                                                                                                                     
  Name:        "IOManager"                                                                                                                           
  IOMode:      1                                                                                                                                     
  OutFileName: "%s"                                                                                                                                  
}                                                                                                                                                    
'''
    #f=open('tmp.cfg','w')                                                                                                                           
    f=tempfile.NamedTemporaryFile(mode='w')
    f.write(cfg % outname)
    f.flush()
    o = larcv.IOManager(f.name)
    o.initialize()
    f.close()
    return o


def get_flow2supera(config_key):

    driver = flow2supera.driver.SuperaDriver()
    if os.path.isfile(config_key):
        driver.ConfigureFromFile(config_key)
    else:
        driver.ConfigureFromFile(flow2supera.config.get_config(config_key))
    
    return driver 

def log_supera_integrity_check(data, driver, log, verbose=False):

    if not log:
        return

    label = driver.Label()
    meta  = driver.Meta()

    # Packet tensor
    pcloud = np.array([[edep.x,edep.y,edep.z,edep.e] for edep in driver._edeps_all])
    voxels = meta.edep2voxelset(driver._edeps_all)
    voxels = np.array([[meta.pos_x(vox.id()),
                        meta.pos_y(vox.id()),
                        meta.pos_z(vox.id()),
                        vox.value()] for vox in voxels.as_vector()
                      ]
                     ) 
    
    cluster_sum = np.sum([p.energy.sum() for p in label.Particles()])
    input_sum  = np.sum([np.sum([edep.e for edep in p.pcloud]) for p in data])
    input_unass= np.sum([edep.e for edep in data.unassociated_edeps])
    energy_sum = label._energies.sum()
    energy_num = label._energies.size()
    pcloud_sum = np.sum(pcloud[:,3])
    pcloud_num = pcloud.shape[0]
    voxels_sum = np.sum(voxels[:,3])
    voxels_num = voxels.shape[0]
    unass_sum  = np.sum([vox.value() for vox in label._unassociated_voxels.as_vector()])
    
    if verbose:
        print('  Raw image    :',voxels_num,'voxels with the sum',voxels_sum)
        print('  Packets      :',pcloud_num,'packets with the sum',pcloud_sum)
        print('  Input cluster:',input_sum)
        print('  Input unass. :',input_unass)
        print('  Label image  :',energy_num,'voxels with the sum',energy_sum)
        print('  Label cluster:',cluster_sum)
        print('  Unassociated :',unass_sum)
        print('  Label image - (Cluster sum + Unassociated)',energy_sum - (cluster_sum + unass_sum))
        print('  Label image - Raw image',energy_sum - voxels_sum)
        print('  Packets - Raw image',pcloud_sum - voxels_sum)

    log['raw_image_sum'].append(voxels_sum)
    log['raw_image_npx'].append(voxels_num)
    log['raw_packet_sum'].append(pcloud_sum)
    log['raw_packet_num'].append(pcloud_num)
    log['in_cluster_sum'].append(input_sum)
    log['in_unass_sum'].append(input_unass)
    log['out_image_sum'].append(energy_sum)
    log['out_image_num'].append(energy_num)
    log['out_cluster_sum'].append(cluster_sum)
    log['out_unass_sum'].append(unass_sum)
    
# Fill SuperaAtomic class and hand off to label-making
def run_supera(out_file='larcv.root',
               in_file='',
               config_key='',
               num_events=-1,
               num_skip=0,
            #    ignore_bad_association=True,
               save_log=None,
               verbose=False):

    is_sim = False
    start_time = time.time()

    writer = get_iomanager(out_file)
  
    driver = get_flow2supera(config_key)
    reader = flow2supera.reader.InputReader(driver.parser_run_config(), in_file,config_key)
    
    if config_key:
        if os.path.isfile(config_key):
            file=config_key
        else:
            file=flow2supera.config.get_config(config_key)
        with open(file,'r') as f:
            cfg=yaml.load(f.read(),Loader=Loader)
            if 'Type' in cfg.keys():
                is_sim=cfg.get('Type')[0]=='sim'
    
    
    

    


    id_vv = ROOT.std.vector("std::vector<unsigned long>")()
    value_vv = ROOT.std.vector("std::vector<float>")()

    id_v=ROOT.std.vector("unsigned long")()
    value_v=ROOT.std.vector("float")()

    if num_events < 0:
        num_events = len(reader)

    print("--- startup {:.2e} seconds ---".format(time.time() - start_time))

    LOG_KEYS  = ['event_id','time_read','time_convert','time_generate', 'time_store', 'time_event']
    LOG_KEYS += ['raw_image_sum','raw_image_npx','raw_packet_sum','raw_packet_num',
    'in_cluster_sum','in_unass_sum','out_image_sum','out_image_num',
    'out_cluster_sum','out_unass_sum']

    logger = dict()
    if save_log:
        for key in LOG_KEYS:
            logger[key]=[]
        driver.log(logger)
    
    print("importing",cfg.get('Type'))

    print("----------------Processing charge events----------------")
    for entry in range(len(reader)):

        if num_skip and entry < num_skip:
            continue

        if num_events <= 0:
            break

        num_events -= 1 

        print(f'Processing Entry {entry}')

        t0 = time.time()
        input_data = reader.GetEntry(entry)

        if input_data.trajectories is None and is_sim:
            print(f'[SuperaDriver] WARNING skipping this entry {entry} as it appears to be "empty" (no truth association found, non-unique event id, etc.)')
            continue

        reader.EventDump(input_data)

        #is_good_event = reader.CheckIntegrity(input_data, ignore_bad_association)
        #if not is_good_event:
        #    print('[ERROR] Entry', entry, 'is not valid; skipping')
        #    continue
        time_read = time.time() - t0
        
        t1 = time.time()
        EventInput = driver.ReadEvent(input_data,is_sim=is_sim)
        time_convert = time.time() - t1

        t2 = time.time()
        driver.GenerateImageMeta(EventInput)
        meta   = larcv_meta(driver.Meta())
        time_generate = time.time() - t2

        # Perform an integrity check
        if save_log:
            log_supera_integrity_check(EventInput, driver, logger, verbose=False)
        t3 = time.time()

        tensor_hits = writer.get_data("sparse3d", "hits")
        if is_sim: 
            driver.Meta().edep2voxelset(driver._edeps_all).fill_std_vectors(id_v, value_v)
        if not is_sim: 
            driver.Meta().edep2voxelset(EventInput.unassociated_edeps).fill_std_vectors(id_v, value_v)

        larcv.as_event_sparse3d(tensor_hits, meta, id_v, value_v)

        if is_sim:
            # Start data store process
            driver.GenerateLabel(EventInput) 
            result = driver.Label()
       
        
            tensor_energy = writer.get_data("sparse3d", "pcluster")
            result.FillTensorEnergy(id_v, value_v)
            larcv.as_event_sparse3d(tensor_energy, meta, id_v, value_v)

            # Check the input image and label image match in the voxel set
            ids_input = np.array([v.id() for v in tensor_energy.as_vector()])
            ids_label = np.array([v.id() for v in tensor_hits.as_vector()])
            assert np.allclose(ids_input,ids_label), '[SuperaDriver] ERROR: the label and input data has different set of voxels'

            tensor_semantic = writer.get_data("sparse3d", "pcluster_semantics")
            result.FillTensorSemantic(id_v, value_v)
            larcv.as_event_sparse3d(tensor_semantic,meta, id_v, value_v)

            cluster_energy = writer.get_data("cluster3d", "pcluster")
            result.FillClustersEnergy(id_vv, value_vv)
            larcv.as_event_cluster3d(cluster_energy, meta, id_vv, value_vv)

            cluster_dedx = writer.get_data("cluster3d", "pcluster_dedx")
            result.FillClustersdEdX(id_vv, value_vv)
            larcv.as_event_cluster3d(cluster_dedx, meta, id_vv, value_vv)

            particle = writer.get_data("particle", "pcluster")
            for p in result._particles:
                if not p.valid:
                    continue
                larp = larcv_particle(p)
                particle.append(larp)

            #Fill mc truth neutrino interactions
            interaction = writer.get_data("neutrino", "mc_truth")
            for ixn in input_data.interactions:
                if isinstance(ixn,np.void):
                    continue
                larn = larcv_neutrino(ixn)
                interaction.append(larn)
            
        #propagating trigger info
        trigger = writer.get_data("trigger", "base")
        trigger.id(int(input_data.event_id))  # fixme: this will need to be different for real data?
        trigger.time_s(int(input_data.t0))
        trigger.time_ns(int(1e9 * (input_data.t0 - trigger.time_s())))   

        # TODO fill the run ID 
        writer.set_id(0, 0, int(input_data.event_id))
        writer.save_entry()

    
   

    
    writer.finalize()

    
    if save_log:
        np.savez('log_flow2supera.npz',**logger)

    end_time = time.time()
    
    print("Total processing time in s: ", end_time-start_time)








