from posixpath import join
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import rc
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple
import argparse
from numpy.core.fromnumeric import size
import scipy.optimize as sco
import networkx as nx
import os,sys,json
import pprint

class PropertySet:
    directory = ""
    properties = dict()
    figsize = tuple()
    filetype = str()
    
    def __init__(self, directory, properties, figsize, filetype):
        self.directory = directory
        self.properties = properties
        self.figsize = figsize
        self.filetype = filetype

def lin(x, a, b):
    return a + b*x

def exp(x, a, b, c):
    return a + b*np.exp(c * x)

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0],tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def plot_time(propertySet):
    print('*  Plotting running time')
    distances = np.array([3,5,7,9,11,13,15,17,21])
    n_qubits = distances**2+(distances-1)**2
    n_qubits_decoder = 2*n_qubits
    n_overlapping_qubits = n_qubits_decoder - 4*distances
    n_messages = 2 * n_overlapping_qubits
    n_checks = n_qubits-1
    
    times = np.array([0.47,4.69,13.48,39.34,75.9,149.0,243.0,343.0,832.0])
    
    p0s = [[-7,3,0.3],[-40,45,0.003],[-40,45,0.005]]
    
    x_data = [distances,n_qubits_decoder,n_messages]
    
    fig,ax = plt.subplots(1,3,sharey=True,figsize=propertySet.figsize,constrained_layout=True)
    ax = ax.ravel()
    for i in range(3):
        ax[i].plot(x_data[i],times,marker='o',linestyle='--')
        popt, pcov = sco.curve_fit(exp,x_data[i],times,p0=p0s[i])
        x_fit = np.linspace(x_data[i][0]-(x_data[i][-1]-x_data[i][0])/10,x_data[i][-1]+(x_data[i][-1]-x_data[i][0])/10,10)
        y_fit = exp(x_fit,*popt)
        ax[i].plot(x_fit,y_fit,marker='',linestyle='-',color='r')
        ax[i].grid()
        
        # print(exp(15,*popt))
        # ax[i].set_ylim(0,160)
    

    
    ax[0].set_xlabel('distance $d$')
    ax[1].set_xlabel('$n_{\mathrm{qubits}}$')
    ax[2].set_xlabel('$n_{\mathrm{messages}}$')
    ax[0].set_ylabel('running time $t/sll$')
    

        
    
    plt.show()

def plot_egs(propertySet):
    print('*  Plotting weight of error guesses, their syndromes and thermodynamic quantities')
    error_guesses = np.load(propertySet.directory+'/errorguesses.npy').astype(np.int32)
    syndromes = np.load(propertySet.directory+'/syndromes.npy').astype(np.int32)
    thermodynamic_quantities = np.load(propertySet.directory+'/tdqs.npy')
    incompatibility_score = np.load(propertySet.directory+'/incompatibility.npy').astype(np.int32)

    weights_e = np.sum(error_guesses,axis=1)
    weights_s = np.sum(syndromes,axis=1)
    weights_e[weights_e<0] = -1
    weights_s[weights_s<0] = -1
    try:
        max_repetitions = int(np.argwhere(weights_s == 0)[0]) + 1
    except:
        max_repetitions = len(weights_e)
    mask = weights_e>=0
    
    avg_energy = thermodynamic_quantities[:,0]
    entropy = thermodynamic_quantities[:,1]
    free_energy = thermodynamic_quantities[:,2]
    
    weights_e = np.array(weights_e[mask])
    weights_s = np.array(weights_s[mask])
    incompatibility_score = np.array(incompatibility_score[mask])
    
    avg_energy = np.array(avg_energy[mask])
    entropy = np.array(entropy[mask])
    free_energy = np.array(free_energy[mask])
    
    max_total_iterations = (propertySet.properties['max_iterations']+1) * propertySet.properties['max_repetitions']
    repetitions = np.zeros(max_total_iterations)
    for i in [propertySet.properties['max_iterations']*j for j in range(1,propertySet.properties['max_repetitions']+1)]:
        repetitions[i] = 1
    
    iterations = np.arange(max_total_iterations)
    iterations = iterations[mask]
    repetitions = repetitions[mask]

    fig,ax = plt.subplots(2,1,sharex=True,figsize=propertySet.figsize,constrained_layout=True)
    ax = ax.ravel()
    ax[0].plot(weights_e,marker='$e$',linestyle=None,color='b')
    ax[0].plot(weights_s,marker='$s$',linestyle=None,color='g')
    ax[0].plot(incompatibility_score,marker='$i$',linestyle=None,color='r')
    ax[1].plot(avg_energy,'bo--',label='$U$')
    ax[1].plot(entropy,'ro--',label='$H$')
    ax[1].plot(free_energy,'go--',label='$F$')
    
    for i,j in enumerate(repetitions):
        tick = repetitions[i]
        if tick > 0 and tick < max_repetitions:
            ax[0].axvline(i-0.5,color='r')
            ax[1].axvline(i-0.5,color='r')

    ax[1].set_xlabel('iteration')
    ax[0].set_ylabel('$|\hat{e}|,|s|,i$')
    ax[1].set_ylabel('Energy')
    ax[1].legend()
    for i in range(2):
        ax[i].grid()
        # ax[i].set_ylim(-0.5)
        # ax[i].set_ylim(-0.5)
    plt.show()

def plot_code(propertySet):
    import matplotlib.animation as animation
    from celluloid import Camera
    import ast
    
    error_guesses_inp = np.load(propertySet.directory+'/errorguesses.npy').astype(np.int32)
    syndromes_inp = np.load(propertySet.directory+'/syndromes.npy').astype(np.int32)
    
    # get initial error and syndrome
    for f in os.listdir(propertySet.directory):
        if f.endswith('.out'):
            OUT_FILE = f
    with open(os.path.join(propertySet.directory,OUT_FILE),'r') as of:
        lines = of.readlines()
    
    for l,line in enumerate(lines):
        if l == 4 and line[0] == 'y':
            start = 0
            for i,c in enumerate(line):
                if c == '=':
                    start = i+2
            y = line[start:-1]
            y = '{' + y + '}'
            initial_error = ast.literal_eval(y)
        if l == 5 and line[0] == 's':
            start = 0
            for i,c in enumerate(line):
                if c == '=':
                    start = i+2
            y = line[start:-1]
            y = '{' + y + '}'

            initial_syndrome_dict = ast.literal_eval(y)
            initial_syndrome = list(initial_syndrome_dict.keys())

    
    n_q = error_guesses_inp.shape[1]
    n_c = syndromes_inp.shape[1]
    
    n_iter = error_guesses_inp.shape[0]
    
    d = int((1 + np.sqrt(2*n_q-1))/2)
    
    
    errors = []
    syndromes = []
    for i in range(n_iter):
        if sum(error_guesses_inp[i])>=0:
            errors += [{k:error_guesses_inp[i,k] for k in range(len(error_guesses_inp[i])) if error_guesses_inp[i,k]!=0}]
            syndromes += [np.argwhere(syndromes_inp[i]==1).flatten()]
            
    fig, ax = plt.subplots(1,1,figsize=(d,d))
    camera = Camera(fig)
    
    if propertySet.filetype == 'pdf':
         ax = surface_code_graph(d,errors[-1],syndromes[-1],initial_error=initial_error,initial_syndrome=initial_syndrome)
         fig.savefig(propertySet.directory+'/last_error_guess.pdf')
    else:
        for i in range(len(syndromes)):
            ax = surface_code_graph(d,errors[i],syndromes[i],initial_error=initial_error,initial_syndrome=initial_syndrome)
            ax.text(0.5,1.02,'iteration {}'.format(i),verticalalignment='center',transform=ax.transAxes)
            camera.snap()
        animation = camera.animate(interval=500,repeat_delay=2000) 
        if propertySet.filetype != '':
            ft = 'gif'
        else:
            ft = propertySet.filetyoe
        animation.save(propertySet.directory+'/decoding_anim.'+ft)
    # fig.savefig('surface_{}.pdf'.format(d))
    return

def surface_code_graph(d,error,syndrome,initial_error=0,initial_syndrome=0):

    n_qubits = d**2+(d-1)**2
    n_checks = n_qubits - 1

    L = 2*d-1
    
    G = nx.grid_2d_graph(L,L)
    pos = {v:np.array([v[0]/L,v[1]/L]) for v in G.nodes()}

    qubits = [n for n in G.nodes() if (n[0]%2==0 and n[1]%2==0) or (n[0]%2==1 and n[1]%2==1)]
    checks = [n for n in G.nodes() if (n[0]%2==0 and n[1]%2==1) or (n[0]%2==1 and n[1]%2==0)]
    X_checks = [n for n in checks if n[1]%2==0]
    Z_checks = [n for n in checks if n[1]%2==1]

    labels_qubits = dict()
    outer_square = 0
    inner_square = d**2
    

    for i,n in enumerate(qubits):
            if n[0]%2==0:
                    labels_qubits[n] = "{}".format(outer_square)
                    outer_square+=1
            else:
                    labels_qubits[n] = "{}".format(inner_square)
                    inner_square+=1
    inv_labels_qubits = {int(v): k for k, v in labels_qubits.items()}

    labels_X_checks = {n:"{}".format(i) for i,n in enumerate(X_checks)}
    labels_Z_checks = {n:"{}".format(i+len(X_checks)) for i,n in enumerate(Z_checks)}

    inv_labels_X_checks = {int(v): k for k, v in labels_X_checks.items()}
    inv_labels_Z_checks = {int(v): k for k, v in labels_Z_checks.items()}

    X_lit_checks = []
    Z_lit_checks = []
    for c in syndrome:
        if c in inv_labels_Z_checks:
            Z_lit_checks += [inv_labels_Z_checks[c]]
        elif c in inv_labels_X_checks:
            X_lit_checks += [inv_labels_X_checks[c]]

    labels=labels_qubits|labels_X_checks|labels_Z_checks
    
    blue_color = '#6666FF'
    red_color = '#FF0000'
    yellow_color = '#FFF700'
    purple_color = '#C900FF'
    white_color = '#FFFFFF'
    green_color = '#40FF00'
    
    error_colors = [red_color,blue_color,yellow_color]
    
    # draw surface code
    nx.draw_networkx_nodes(G, pos=pos,nodelist=qubits, node_shape='o',edgecolors=purple_color,node_color=white_color,linewidths=2.0)
    nx.draw_networkx_nodes(G, pos=pos,nodelist=Z_checks, node_shape='s',edgecolors=blue_color,node_color=white_color,linewidths=2.0)
    nx.draw_networkx_nodes(G, pos=pos,nodelist=X_checks, node_shape='s',edgecolors=red_color,node_color=white_color,linewidths=2.0)
    
    # draw initial error (if given)
    if initial_error != 0:
        for e in [1,2,3]:
            error_qubits = [inv_labels_qubits[q] for q in initial_error.keys() if initial_error[q]==e]
            nx.draw_networkx_nodes(G, pos=pos,nodelist=error_qubits, node_shape='o',edgecolors=error_colors[e-1],node_color=white_color,alpha=0.5,node_size=500,linewidths=3.0)
    
    # draw initial syndrome (if given)
    if initial_syndrome !=0:
        X_initial_checks = []
        Z_initial_checks = []
        for c in initial_syndrome:
            if c in inv_labels_Z_checks:
                Z_initial_checks += [inv_labels_Z_checks[c]]
            elif c in inv_labels_X_checks:
                X_initial_checks += [inv_labels_X_checks[c]]
        nx.draw_networkx_nodes(G, pos=pos,nodelist=Z_initial_checks, node_shape='s',edgecolors=green_color,node_color=white_color,alpha=0.5,node_size=500,linewidths=3.0)
        nx.draw_networkx_nodes(G, pos=pos,nodelist=X_initial_checks, node_shape='s',edgecolors=green_color,node_color=white_color,alpha=0.5,node_size=500,linewidths=3.0)
    
    
    # draw errors
    for e in [1,2,3]:
        error_qubits = [inv_labels_qubits[q] for q in error.keys() if error[q]==e]
        nx.draw_networkx_nodes(G, pos=pos,nodelist=error_qubits, node_shape='o',edgecolors=purple_color,node_color=error_colors[e-1],linewidths=2.0)
    
    
    # draw lit checks
    nx.draw_networkx_nodes(G, pos=pos,nodelist=Z_lit_checks, node_shape='s',edgecolors=blue_color,node_color=green_color,linewidths=2.0)
    nx.draw_networkx_nodes(G, pos=pos,nodelist=X_lit_checks, node_shape='s',edgecolors=red_color,node_color=green_color,linewidths=2.0)
    

    # draw labels and edges
    nx.draw_networkx_labels(G,pos=pos,labels=labels)
    nx.draw_networkx_edges(G,pos=pos)
    return plt.gca()
    

def main():
    np.set_printoptions(linewidth=sys.maxsize)
    pp = pprint.PrettyPrinter()

    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 11
    rcParams['text.usetex'] = True

    parser = argparse.ArgumentParser(description='Tools to plot results of simulations')
    parser.add_argument('-t','--type',default='egs',help='type of plots to show')
    parser.add_argument('-c','--code_type',default='surface',help='codetype')
    parser.add_argument('-d','--dir',default='tmp',help='output directory')
    parser.add_argument('-i','--n_iter',default=0,help='number of iterations to plot')
    parser.add_argument('-f','--file_type',default='png',help='filetype')
    parser.add_argument('-p','--plot',default='False',help='whether to plot')
    parser.add_argument('-ps','--plot_size',default='10x10',help='plot size in cm')
    parser.add_argument('-po','--plot_order',default='square',help='plot order, square, linear')


    args = parser.parse_args()

    type = args.type
    code_type = args.code_type
    directory = args.dir
    if directory == 'tmp':
        directory = '/Users/josiasold/Documents/Studium/Archiv/ma/ma-paper/ma-paper_numerics/gbp/sim/output/tmp'
    plot = args.plot
    n_iter = int(args.n_iter)
    file_type = args.file_type
    plot_size = args.plot_size

    if args.plot_size.split('x')[0] == 'fw':
        plot_size_w = 11.80737
        plot_size_h = float(args.plot_size.split('x')[1])
    elif args.plot_size.split('x')[0] == 'hw':
        plot_size_w = 11.80737 / 2.0
        plot_size_h = float(args.plot_size.split('x')[1])
    else:
        plot_size_w = float(args.plot_size.split('x')[0])
        plot_size_h = float(args.plot_size.split('x')[1])
    fs = cm2inch(plot_size_w,plot_size_h)
    plot_order = args.plot_order
    
    properties = dict()
    json_files = [f for f in os.listdir(directory) if f.endswith('.json')]
    
    with open(os.path.join(directory,json_files[0])) as property_file:
        properties = json.load(property_file)

    print('** Data directory: {} **'.format(directory))
    print('*  {}:'.format(json_files[0]))
    pp.pprint(properties)
    
    propertySet = PropertySet(directory,properties,fs,file_type)
    
    if type == 'egs':
        plot_egs(propertySet)
    elif type == 'code':
        plot_code(propertySet)
    elif type == 'time':
        plot_time(propertySet)


if __name__ == "__main__":
    main()