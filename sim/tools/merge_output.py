import numpy as np
import os,sys,json,shutil
import argparse

def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    from https://stackoverflow.com/questions/3041986/apt-command-line-interface-like-yes-no-input/3041990
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")

def merge_data(dir1,dir2):
        dir_1 = os.listdir(dir1)
        dir_2 = os.listdir(dir2)
        input_file_1, output_file_1 = [x for x in dir_1 if x.endswith('json') or x.endswith('out')]
        input_file_2, output_file_2 = [x for x in dir_2 if x.endswith('json') or x.endswith('out')]
        print(input_file_1, output_file_1)
        print(input_file_2, output_file_2)
        # file_1 = 
        # os.path.join(dir1,)

        with open(os.path.join(dir1,input_file_1)) as FILE:
                input_dict_1 = json.load(FILE)
        
        with open(os.path.join(dir2,input_file_2)) as FILE:
                input_dict_2 = json.load(FILE)

        assert input_dict_1 == input_dict_2, "Runs are not equivalent, {}".format({k: (input_dict_1[k],input_dict_2[k]) for k in input_dict_1 if k in input_dict_2 and input_dict_1[k] != input_dict_2[k]})
        
        output_dict = input_dict_1
        output_dict['n_errorsamples'] += input_dict_2['n_errorsamples']
        
        FILE_1 = os.path.join(dir1,output_file_1)
        output_data_1 = np.genfromtxt(FILE_1,skip_footer=2,skip_header=3)
        total_time_1 = np.genfromtxt(FILE_1,skip_header=15,max_rows=1)[4]
        
        FILE_2 = os.path.join(dir2,output_file_2)
        output_data_2 = np.genfromtxt(FILE_2,skip_footer=2,skip_header=3)
        total_time_2 = np.genfromtxt(FILE_2,skip_header=15,max_rows=1)[4]
        
        print(output_data_1)
        print(total_time_1)
        print(output_data_2)
        print(total_time_2)
        
        sum_data = output_data_1+output_data_2
        for column in [0,-1,-2,-3]:
                sum_data[:,column]/=2
        
        n_data_rows,n_data_columns = sum_data.shape
        
        OUTDIR = dir1
        # os.makedirs(OUTDIR,exist_ok=True)
        
        with open(os.path.join(OUTDIR,input_file_1),'w+') as OUTPUT_DICT:
                json.dump(output_dict,OUTPUT_DICT,separators=(',\n', ': '))
        
        string_around = []
        with open(FILE_1) as INPUT_OUT:
                for line_number,line in enumerate(INPUT_OUT):
                        if line_number in [0,1,2,n_data_rows+3,n_data_rows+4]:
                                string_around += [line]
        
        print(string_around)
        
        with open(os.path.join(OUTDIR,output_file_1),'w+') as OUTPUT_OUT:
                data_row = 0
                for line_number in range(n_data_rows+5):
                        if line_number <= 2:
                                OUTPUT_OUT.write(string_around[line_number])
                        elif line_number > 2 and data_row < n_data_rows:
                                OUTPUT_OUT.write(np.array2string(sum_data[data_row],separator='\t',formatter={'float_kind':lambda x: "%.5g" % x})[1:-1])
                                OUTPUT_OUT.write('\n')
                                data_row+=1
                        elif line_number > n_data_rows+2 and line_number < n_data_rows+4:
                                OUTPUT_OUT.write(string_around[-2])
                        else:
                                OUTPUT_OUT.write(string_around[-1][:21])
                                OUTPUT_OUT.write('{:.0f}'.format(total_time_1+total_time_2))
                                OUTPUT_OUT.write(' s')
        if query_yes_no("Delete directory {}".format(dir2)):
                shutil.rmtree(dir2)

def main():
        parser = argparse.ArgumentParser(description='Tools to merge outputs of simulations')
        parser.add_argument('-d1','--dir_1',help='input directory 1')
        parser.add_argument('-d2','--dir_2',help='input directory 2')
        
        args = parser.parse_args()
        
        np.set_printoptions(linewidth=sys.maxsize)

        base_dir = '$HOME/Dokumente/gbp_paper/gbp/sim/output/surface'
        mergers = [('1126_132028_surface_11','1126_183528_surface_11'),('1127_025547_surface_11','1127_122412_surface_11'),('1127_170912_surface_11','1127_233311_surface_11'),('1128_082557_surface_11','1128_124953_surface_11')]

        for merge in mergers:
                dir_1 = os.path.join(base_dir,merge[0])
                dir_2 = os.path.join(base_dir,merge[1])

                merge_data(dir_1,dir_2)
        
if __name__ == '__main__':
        main()