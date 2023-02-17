import sys, os
import h5py
import numpy as np
import time
import larnd2supera
import argparse
import ROOT
from edep2supera.utils import get_iomanager, larcv_meta, larcv_particle
from LarpixParser import event_parser as EventParser
from larcv import larcv
import pandas as pd

def get_larnd2supera(config_key):

    driver = larnd2supera.driver.SuperaDriver()
    if os.path.isfile(config_key):
        driver.ConfigureFromFile(config_key)
    else:
        driver.ConfigureFromFile(larnd2supera.config.get_config(config_key))
    
    return driver 


# Fill SuperaAtomic class and hand off to label-making
def run_supera(out_file='larcv.root',
               in_file='',
               config_key='',
               num_events=-1,
               num_skip=0,
               save_log=None):

    start_time = time.time()

    reader = larnd2supera.reader.InputReader(in_file)
    writer = get_iomanager(out_file)
    driver = get_larnd2supera(config_key)

    id_vv=ROOT.std.vector("std::vector<unsigned long>")()
    value_vv=ROOT.std.vector("std::vector<float>")()

    id_v=ROOT.std.vector("unsigned long")()
    value_v=ROOT.std.vector("float")()

    if num_events < 0:
        num_events = len(reader)

    print("--- startup {:.2e} seconds ---".format(time.time() - start_time))

    LOG_KEYS=['event_id','time_read','time_convert','time_generate', 'time_store', 'time_event']
    logger = dict()
    if save_log:
        for key in LOG_KEYS:
            logger[key]=[]
        driver.log(logger)
        
    for entry in range(len(reader)):

        if num_skip and entry < num_skip:
            continue

        if num_events <= 0:
            break

        num_events -= 1 

        t0 = time.time()
        input_data = reader.GetEntry(entry)
        time_read = time.time() - t0

        print('Entry {} Event ID {}'.format(entry,input_data.event_id))
        
        t1 = time.time()
        EventInput = driver.ReadEvent(input_data)
        time_convert = time.time() - t1

        # TODO Seems to run, but how to check it's really working?
        # TODO Should this be supera_driver or driver?
        t2 = time.time()
        driver.GenerateImageMeta(EventInput)
        driver.GenerateLabel(EventInput) 
        time_generate = time.time() - t2

        t3 = time.time()
        value_v  = ROOT.std.vector("float")()
        value_vv = ROOT.std.vector("std::vector<float>")()

        result = driver.Label()
        meta   = larcv_meta(driver.Meta())
        
        tensor_energy = writer.get_data("sparse3d","pcluster")
        result.FillTensorEnergy(id_v,value_v)
        larcv.as_event_sparse3d(tensor_energy,meta,id_v,value_v)

        tensor_semantic = writer.get_data("sparse3d","pcluster_semantics")
        result.FillTensorSemantic(id_v,value_v)
        larcv.as_event_sparse3d(tensor_semantic,meta,id_v,value_v)

        cluster_energy = writer.get_data("cluster3d","pcluster")
        result.FillClustersEnergy(id_vv,value_vv)
        larcv.as_event_cluster3d(cluster_energy,meta,id_vv,value_vv)

        cluster_dedx = writer.get_data("cluster3d","pcluster_dedx")
        result.FillClustersdEdX(id_vv,value_vv)
        larcv.as_event_cluster3d(cluster_dedx,meta,id_vv,value_vv)

        #
        # Unassociated EDeps: both image and semantics should take care of this
        #
        edep_unass = driver._edeps_unassociated
        
        particle = writer.get_data("particle","pcluster")
        for p in result._particles:
            if not p.valid:
                continue
            larp = larcv_particle(p)
            particle.append(larp)

        # TODO fill the run ID 
        writer.set_id(0,0,int(input_data.event_id))
        writer.save_entry()
        time_store = time.time() - t3

        time_event = time.time() - t0
        print("--- running driver  {:.2e} seconds ---".format(time_event))

        if save_log:
            logger['event_id'].append(input_data.event_id)
            logger['time_read'    ].append(time_read)
            logger['time_convert' ].append(time_convert)
            logger['time_generate'].append(time_generate)
            logger['time_store'   ].append(time_store)
            logger['time_event'   ].append(time_event)

    writer.finalize()

    # store supera log dictionary
    if save_log:
        np.savez('log_larnd2supera.npz',**logger)

    print("done")   








