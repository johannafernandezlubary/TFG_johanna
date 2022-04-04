print("porfa")
import sys
import numpy as np
import bluesky as bs
from bluesky import traffic as tr
from bluesky import settings
from bluesky.traffic.route import Route
from bluesky.navdatabase import Navdatabase
from bluesky.simulation import Simulation
from bluesky.traffic.performance.perfbase import PerfBase
import matplotlib.pyplot as plt
from bluesky.tools import geo
import random
from sacred import Experiment
import networkx as nx
from bluesky.traffic.asas.mvp import MVP as mvp
print("passa los imports")
      
ex = Experiment("test-experiment")

print ("passa ex")

def check_boundaries(traf, center, radius):
    """
    Check if any aircraft is out of the scenario bounds. It deletes it if so.
    """
    radius = radius * 1852 # From nm to meters
    id_to_delete = []
    for i in range(traf.ntraf):
        if geo.latlondist(traf.lat[i], traf.lon[i] , center[0], center[1]) > radius:
            id_to_delete.append(traf.id[i])

    if id_to_delete:
        for idx in id_to_delete:
            traf.delete(bs.traf.id.index(idx))
            

def create_sources(center, radius, n_sources):
    """
    The sources will create a polygon centered in the simulation. The function returns
    a list with each source coordinates
    """
    sources_positions = []

    dist_to_center = radius

    for i in range(n_sources):
        alpha = 360/n_sources * i
        lat, lon = geo.qdrpos(center[0], center[1], alpha, dist_to_center)
        sources_positions.append([lat, lon])

    return sources_positions

def init_at(radius, center, n_ac):

    earth_radius = geo.rwgs84(center[0])/1852
    ac_list=[]
    for _ in range(n_ac):

        random_angle = random.random() * 360
        random_distance = random.random()

        random_lat, random_lon = geo.qdrpos(center[0], center[1], random_angle, random_distance)


        angle = geo.qdrdist(random_lat, random_lon, center[0], center[1])[0]
        limit_angle = 60.0 # 60deg when the sources are in the np

        acid = str(random.getrandbits(32))
        heading = random.uniform(angle - limit_angle, angle + limit_angle)
        speed=random.uniform(38.8768898, 42.7645788 )
        alt= 49.2125984 #15 m to ft
        actype="M200"
        
        #bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon,  achdg=heading, acalt=alt, acspd=speed )
        ac_list.append([acid, actype, random_lat, random_lon, heading, alt, speed])
        
        return ac_list
    
def spawn_ac(sources_position, radius, center, number_of_aircrafts):

    for _ in range(number_of_aircrafts):
        source = random.choice(sources_position)
        angle = geo.qdrdist(source[0], source[1], center[0], center[1])[0]
        limit_angle = 60.0 * 180 / np.pi

        acid = str(random.getrandbits(32))
        heading = random.uniform(angle - limit_angle, angle + limit_angle)
        
        speed=random.uniform(38.8768898, 42.7645788 )
        alt= 49.2125984 #15 m to ft

        bs.traf.cre(acid, actype="M200", aclat=source[0], aclon=source[1], achdg=heading, acalt=alt, acspd=speed)


def plot_at(center, radius, sources_position):

    alpha = np.linspace(0, 360, 500)


    coords = [geo.qdrpos(center[0], center[1], angle, radius) for angle in alpha]

    x = bs.traf.lat
    y = bs.traf.lon
    vx = bs.traf.gseast
    vy = bs.traf.gsnorth

    plt.figure(figsize=(8, 8))
    plt.scatter(center[0], center[1], 2)
    plt.scatter([coord[0] for coord in coords], [coord[1] for coord in coords], 1)
    plt.scatter([coord[0] for coord in sources_position], [coord[1] for coord in sources_position])
    plt.quiver(x, y, vy, vx)
    plt.title("Scenario initial configuration")
    plt.show()

def comp_conf(graph):
        """
        Return a list with all compound conficts
        """
        if graph.number_of_edges != 0:
            comp_confs = [comp_conf for comp_conf in nx.connected_components(graph) if len(comp_conf) > 2] 
            return comp_confs

        return []    
    
    
    
def simstep():
    bs.sim.step()
    bs.net.step()

def change_ffmode(mode=True):
        bs.sim.ffmode = mode
        bs.sim.dtmult = 100.0



def at_to_graph():
    graph = nx.Graph()
    graph.add_nodes_from(bs.traf.id)
    for pair in bs.traf.cd.confpairs_unique:
        pair=list(pair)
        graph.add_edge(pair[0], pair[1])
     
    #print("entra en at to graph")
    
    return graph
  

    
def append_variables(variables, graph):
    """
    During a time window when conflicts are present stores the values of the variables
    in the dictionary so maximum values can be computed later
    """
    for pair in bs.traf.cd.lospairs_unique:
        pair=list(pair)
        pair = {pair[0], pair[1]}
        if (not pair in variables["los"]):
            variables["los"].append(pair)
            #print(variables["los"])
            
    for pair in bs.traf.cd.confpairs_unique:
        pair = list(pair)
        pair = {pair[0], pair[1]}
        if (not pair in variables["confs"]):
            variables["confs"].append(pair)
            
    comp_confs = comp_conf(graph)
    for conf in comp_confs:
        if not conf in variables["comp_confs"]:
            variables["comp_confs"].append(conf)            
            
   
def log_noCR(variables, _run):
    _run.log_scalar("number_los_noCR", len(variables["los"]))

    _run.log_scalar("number_conflicts_noCR", len(variables["confs"]))
    
    _run.log_scalar("number_comp_conf_noCR",len(variables["comp_confs"]))

    
    if variables["comp_confs"]:
        conf_size = max([len(comp_conf) for comp_conf in variables["comp_confs"]])
    else:
        conf_size = 0

    _run.log_scalar("conf_size_noCR", conf_size)
    
def log_withCR(variables, _run):
    _run.log_scalar("number_los_withCR", len(variables["los"]))

    _run.log_scalar("number_conflicts_withCR", len(variables["confs"]))
    

    _run.log_scalar("number_comp_conf_withCR",len(variables["comp_confs"]))

    
    if variables["comp_confs"]:
        conf_size = max([len(comp_conf) for comp_conf in variables["comp_confs"]])
    else:
        conf_size = 0

    _run.log_scalar("conf_size_withCR", conf_size)
    
    
@ex.config
def cfg():
    
    center = (47, 9)
    radius = 1.34989201 #nautical miles change it
    n_ac = 50
    sim_time =300 #seconds = 15min
    n_sources = 10
    n_runs = 3

@ex.automain
def complexity_simulation(_run, center, radius, n_ac, sim_time, n_sources, n_runs):
    

    bs.init('sim-detached')
    
      
    #bs.traf.cd.rpz_def =0.0890928726 #165meters into nauticalmiles/  0.129589633natuical miles are 240meters
    #bs.traf.cd.dtlookahead =15.0  #seconds  TELL RALVI TO CHANGE DT
    #bs.settings.asas_dtlookahead=15
    
    #bs.traf.cd.setmethod("ON")
    #bs.traf.cr.setmethod('MVP')


    print("inicia simu")
    
    
    for run in range(n_runs):
        
        control=True
        sources_position = create_sources(center, radius, n_sources)
        ac_list=init_at(radius, center, n_ac) #here calls the function to create the aircrafts
        
        if control==True:
            print("ini simu sin CR")
            _run.log_scalar("randomvalue_noCR", 89)
            bs.traf.cd.setmethod("ON")
            #bs.traf.cr.setmethod('MVP')

            bs.traf.cd.rpz_def =0.0890928726 #165meters into nauticalmiles/  0.129589633natuical miles are 240meters
            bs.traf.cd.dtlookahead =15.0

            #bs.sim.simdt = 1
           # bs.sim.simt = 0

            t_max = sim_time #15 mins

            ntraf = bs.traf.ntraf
            n_steps = int(t_max//bs.sim.simdt + 1)
            print(n_steps)
            t = np.linspace(0, t_max, n_steps)    
            
            for i in range(len(ac_list)):
                #bs.traf.cre(int(ac_list[i][0]),str(ac_list[i][1]),int(ac_list[i][2]),int(ac_list[i][3]),int(ac_list[i][4]),int(ac_list[i][5]),int(ac_list[i][6]))
                #bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon,  achdg=heading, acalt=alt, acspd=speed )
                acid = ac_list[i][0]
                actype = ac_list[i][1]
                aclat = ac_list[i][2]
                aclon = ac_list[i][3]
                achdg = ac_list[i][4]
                acalt = ac_list[i][5]
                acspd = ac_list[i][6]
                
                bs.traf.cre(acid, actype, aclat, aclon,  achdg, acalt, acspd)


            variables = {"confs": [], "los":[], "comp_confs": [], "rndm":[]}
            """ Main loop """
            print("loop noCR")
            
            change_ffmode()
            
            print("ff")
            for i in range(n_steps):
                #print(bs.traf.cd.lospairs_unique)


                """ Check if the acs are out of bounds and delete them if so """
                check_boundaries(bs.traf, center, radius)

                """ Spawning aircrafts in the sources """
                if bs.traf.ntraf < n_ac:
                    spawn_ac(sources_position, radius, center, number_of_aircrafts = n_ac - bs.traf.ntraf)

                graph = at_to_graph()
                append_variables(variables, graph)
                simstep()
                
            log_noCR(variables,_run)
            control=False
        
        if control==False:
            
            print("ini simu CR")
            _run.log_scalar("randomvalue_withCR", 7777)
            bs.traf.cd.setmethod("ON")
            bs.traf.cr.setmethod('MVP')

            bs.traf.cd.rpz_def =0.0890928726 #165meters into nauticalmiles/  0.129589633natuical miles are 240meters
            bs.traf.cd.dtlookahead =15.0

            #bs.sim.simdt = 1
           # bs.sim.simt = 0

            t_max = sim_time #15 mins

            ntraf = bs.traf.ntraf
            n_steps = int(t_max//bs.sim.simdt + 1)
            print(n_steps)
            t = np.linspace(0, t_max, n_steps)    
            
            for i in range(len(ac_list)):
                #bs.traf.cre(int(ac_list[i][0]),str(ac_list[i][1]),int(ac_list[i][2]),int(ac_list[i][3]),int(ac_list[i][4]),int(ac_list[i][5]),int(ac_list[i][6]))
                #bs.traf.cre(acid, actype="M200", aclat=random_lat, aclon=random_lon,  achdg=heading, acalt=alt, acspd=speed )
                acid = ac_list[i][0]
                actype = ac_list[i][1]
                aclat = ac_list[i][2]
                aclon = ac_list[i][3]
                achdg = ac_list[i][4]
                acalt = ac_list[i][5]
                acspd = ac_list[i][6]
                
                bs.traf.cre(acid, actype, aclat, aclon,  achdg, acalt, acspd)
                
            variables = {"confs": [], "los":[], "comp_confs": [], "rndm":[]}
            """ Main loop """
            print("loop noCR")
            
            change_ffmode()
            
            print("ff")
            for i in range(n_steps):
                #print(bs.traf.cd.lospairs_unique)


                """ Check if the acs are out of bounds and delete them if so """
                check_boundaries(bs.traf, center, radius)

                """ Spawning aircrafts in the sources """
                if bs.traf.ntraf < n_ac:
                    spawn_ac(sources_position, radius, center, number_of_aircrafts = n_ac - bs.traf.ntraf)

                graph = at_to_graph()
                append_variables(variables, graph)
                simstep()
                
            log_withCR(variables,_run)
            control=True
    bs.sim.reset()   