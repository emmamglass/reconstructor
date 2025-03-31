import os

def reconstruct(input_file, file_type = 1, media=[], org = 'default', min_frac = 0.01, max_frac = 0.5, gram='none', out = 'default', name = 'default', cpu = 1, gapfill = 'yes', test = 'no'):
    print('Generating reconstruction....')
    cmd_line = 'python -m reconstructor --input_file %s --file_type %s --media %s --org %s --min_frac %s --max_frac %s --gram %s --out %s --name %s --cpu %s --gapfill %s --test %s'%(input_file,file_type,media,org,min_frac,max_frac,gram,out,name,cpu,gapfill,test)
    print('Reconstruction generating with the following command line:')
    print(cmd_line)
    os.system(cmd_line)

if __name__ == '__main__':
    reconstruct('218496.4.fa', file_type = 1, gram = 'negative')