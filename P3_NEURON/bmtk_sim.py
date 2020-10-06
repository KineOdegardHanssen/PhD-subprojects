

def build_sim():
    from bmtk.builder.networks import NetworkBuilder
    
    # My all active-model (does not work):
    #Model ID 496497595
    #Cell ID 487667205
    
    # Other perisomatic model (available on Allen Brain Institute - CellTypes):
    #Model ID 491623973
    #Cell  ID 490387590
    
    print('BMTK import success')
    net = NetworkBuilder('mcortex')
    print('Network builder initiated')
    # Testing other models: # Surprise, surprise, this does not work...
    net.add_nodes(cell_name='Pvalb_490387590_m',
                  model_type='biophysical',
                  model_template='ctdb:Biophys1.hoc',
                  model_processing='aibs_perisomatic',
                  dynamics_params='491623973_fit.json',
                  morphology='Pvalb_490387590_m.swc')
    ''' # Standard
    net.add_nodes(cell_name='Scnn1a_473845048',
                  potental='exc',
                  model_type='biophysical',
                  model_template='ctdb:Biophys1.hoc',
                  model_processing='aibs_perisomatic',
                  dynamics_params='472363762_fit.json',
                  morphology='Scnn1a_473845048_m.swc')
    '''
    print('Node added')

    net.build()
    net.save_nodes(output_dir='network')
    for node in net.nodes():
        print(node)
    print('Node printed')


    from bmtk.utils.sim_setup import build_env_bionet
    print('Setting environment')

    build_env_bionet(base_dir='sim_ch01',      # Where to save the scripts and config files
                     network_dir='network',    # Location of directory containing network files
                     tstop=1200.0, dt=0.1,     # Run a simulation for 2000 ms at 0.1 ms intervals
                     report_vars=['v'], # Tells simulator we want to record membrane potential and calcium traces
                     current_clamp={           # Creates a step current from 500.ms to 1500.0 ms
                         'amp': 0.61,   # 0.12#0.610
                         'delay': 100.0, # 100, #500
                         'duration': 1000.0
                     },
                     include_examples=True,    # Copies components files
                     compile_mechanisms=True   # Will try to compile NEURON mechanisms
                    )
    print('Build done')


def run_sim():
    from bmtk.simulator import bionet


    conf = bionet.Config.from_json('sim_ch01/simulation_config.json')
    conf.build_env()
    net = bionet.BioNetwork.from_config(conf)
    sim = bionet.BioSimulator.from_config(conf, network=net)
    sim.run()

    from bmtk.analyzer.spike_trains import to_dataframe
    to_dataframe(config_file='sim_ch01/simulation_config.json')

    from bmtk.analyzer.cell_vars import plot_report

    plot_report(config_file='sim_ch01/simulation_config.json')


if __name__ == '__main__':
    #build_sim()
    run_sim()