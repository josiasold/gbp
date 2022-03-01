import numpy as np
import os,time,json,sys
import argparse
import itertools as it
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pprint
from shutil import get_terminal_size
import scipy.optimize as sco

def lin(x, a, b):
    return a + b*x

def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0],tuple):
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)


def main():
    parser = argparse.ArgumentParser(description='Tools to plot results of simulations')
    parser.add_argument('-d','--dir',default='latest',help='input directory')
    parser.add_argument('-s','--save_fig',default='True',help='whether figure(s) shall be saved')
    parser.add_argument('-x','--xaxis',default='p_error',help='what to put on the x axis')
    parser.add_argument('-y','--yaxis',default='ber',help='what to put on the y axis')
    parser.add_argument('-p','--parameter',default='code',help='which parameter to use for different curves')
    parser.add_argument('-t','--tex',default='True',help='whether TeX should be used to render figures')
    parser.add_argument('-ps','--plot_size',default='20x20',help='Size of plot in cm: widthxheight e.g. 20x10 for widht 20cm height 10cm')
    parser.add_argument('-st','--statistics',default='1',help='whether to consider statistics, and how many runs per point')
    parser.add_argument('-f','--fit',default='False',help='whether to fit a linear function to the data points')
    
    args = parser.parse_args()

    pp = pprint.PrettyPrinter(indent=4,width=get_terminal_size()[0],depth=4)
    np.set_printoptions(suppress=True,linewidth=sys.maxsize,threshold=sys.maxsize)

    if args.tex == 'True':
        rcParams['text.usetex'] = True
    else:
        rcParams['text.usetex'] = False
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 11
    
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

    
    DIRECTORY ="output/"+args.dir
    colors = {}
    
    plot_x = args.xaxis
    plot_y = args.yaxis
    parameter = args.parameter
    save_fig = args.save_fig
    n_statistics = int(args.statistics)
    fit = args.fit
    
    data = dict()

    for path,dirs,files in os.walk(DIRECTORY):
        if len(files) > 1:
            for file in files:
                if file.endswith('.out'):
                    code = file.rsplit('.')[0]
                    outfile = os.path.join(path,file)
                if file.endswith('.json'):
                    jsonfile = os.path.join(path,file)
            if code not in data:
                data[code] = dict()
                for key in ["p_error","ch_sc","sci","sce","ler","fail","repeat","iter","ber"]:
                    data[code][key] = []
            
            n_errorsamples = sum(np.genfromtxt(outfile,delimiter='\t',skip_header=3,max_rows=1)[1:6])
            input_data = np.genfromtxt(outfile,delimiter='\t',skip_header=3,skip_footer=2)

            data[code]["p_error"] += [input_data[:,0]]
            data[code]["ch_sc"] += [input_data[:,1]/n_errorsamples]
            data[code]["sci"] += [input_data[:,2]/n_errorsamples]
            data[code]["sce"] += [input_data[:,3]/n_errorsamples]
            data[code]["ler"] += [input_data[:,4]/n_errorsamples]
            data[code]["fail"] += [input_data[:,5]/n_errorsamples]
            data[code]["repeat"] += [input_data[:,6]] 
            data[code]["iter"] += [input_data[:,7]]
            data[code]["ber"] += [input_data[:,-1]]


    fig,ax = plt.subplots(1,1,figsize=fs,constrained_layout=True)

    for code in data:
        for key in data[code]:
            data[code][key] = np.array(data[code][key])

    for code in data:
        print(data[code])
        try:
            x_data = np.mean(np.array(data[code][plot_x]),axis=0)
            y_data = np.mean(np.array(data[code][plot_y]),axis=0)
            y_std = np.std(np.array(data[code][plot_y]),axis=0,ddof=1) #/ np.sqrt(np.array(data[code][plot_y]).shape[0])
            
            linestyle = '--'
            if fit == 'lin':
                linestyle = ''
                popt, pcov = sco.curve_fit(lin,x_data,y_data,sigma=y_std)
                x_fit = np.linspace(x_data[0]-(x_data[-1]-x_data[0])/10,x_data[-1]+(x_data[-1]-x_data[0])/10,3)
                y_fit = lin(x_fit,*popt)
            
            color = next(ax._get_lines.prop_cycler)['color']
            ax.errorbar(x_data,y_data,yerr=y_std,marker='x', ls = linestyle,label=code.rsplit('_')[-1],color=color)
            ax.plot(x_fit,y_fit,color=color)
        except:
            print("Fail")
            continue

        
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left', ncol=len(data.keys()), mode="expand", borderaxespad=0.,handlelength=1.0,handletextpad = 0.01,columnspacing = 1.0, title = 'distance')
    ax.grid()
    # ax.set_ybound(0)
    if plot_x == 'p_error':
        ax.set_xlabel('$p_{\mathrm{error}}$')
    if plot_y == 'ch_sc':
        ax.set_ylabel('$p_{\mathrm{channel\:success}}$')
    elif plot_y == 'ber':
        ax.set_ylabel('$p_{\mathrm{block\:error}}$')
    elif plot_y == 'ler':
        ax.set_ylabel('$p_{\mathrm{logical\:error}}$')
    elif plot_y == 'sci':
        ax.set_ylabel('$p_{\mathrm{success\:identical}}$')
    elif plot_y == 'sce':
        ax.set_ylabel('$p_{\mathrm{success\:equivalent}}$')
    elif plot_y == 'repeat':
        ax.set_ylabel('average $\#$ of repetitions')
    elif plot_y == 'iter':
        ax.set_ylabel('average $\#$ of iterations')
    elif plot_y == 'fail':
        ax.set_ylabel('$p_{\mathrm{failure}}$')
    

    if save_fig == 'True':
        DIR = os.path.join(DIRECTORY,'plot_{}_{}_{}.pdf'.format(plot_x,plot_y,parameter))
        fig.savefig(DIR)
    plt.show()
    
if __name__ == '__main__':
    main()