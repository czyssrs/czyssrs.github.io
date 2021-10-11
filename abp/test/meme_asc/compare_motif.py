#
# Fri Mar 17 00:49:34 PDT 2017
# Author: Honglei Liu <liuhonglei@gmail.com>
#
# Description: merge or compare between pwms
#


import sys
import os
import math
import parse_motif

from yattag import Doc

"'using KL-divergence to find similar motifs and merge them together'"
#bg=[]           # background model
kl_bound=1.0      # KL threshold to merge two motifs
beta=0.00001    # remove zeros
minL=4      # segement length
e_thres=0.01  # threshold for e-value
min_occ=10  # minmum number of occurrences

def calc_kl(m1,m2):
    '''m2 shifts, using KL-divergence to compare'''

    #if(x1['nsites']>x2['nsites']):
    #    m1,m2=x1,x2
    #else:
    #    m1,m2=x2,x1

    if m1['alength']!=m2['alength']:
        print "pwms have different Alphabet set!"
        return

    pwm1=m1['pwm_p']
    pwm2=m2['pwm_p']
    lag=m2['w']-m1['w']
    p_value=sys.float_info.max
    shift=0

    for i in range(0-(m2['w']-minL),m1['w']-minL+1):     #at least have 6 overlaps
        t_p=0
        l_1=max(0,i)
        r_1=min(m1['w'],m1['w']+min(0,lag+i))
        l_2=max(0,0-i)
        r_2=min(m2['w'],m2['w']-max(0,lag+i))

        for j in range(r_1-l_1):
            for k in range(m1['alength']):
                # make sure the probability is not zero
                tp1=pwm1[l_1+j][k]
                if tp1==0:
                    tp1+=beta
                tp2=pwm2[l_2+j][k]
                if tp2==0:
                    tp2+=beta
                # calculate KL-Divergence
                t_p+=tp1*(math.log(tp1)-math.log(tp2))
        t_p=t_p/(r_1-l_1)
        if t_p<p_value:
            p_value=t_p
            shift=i

    return p_value,shift


def compare(par1, logo_dir1, par2, logo_dir2):

    print "[1] motifs number: "+str(len(par1.data))
    print "[2] motifs number: "+str(len(par2.data))
    common=0
    logos, meta= [], []
    for i in range(len(par1.data)):
        min_kl = float('inf')
        j = 0
        shift = 0
        for k in range(len(par2.data)):
            kl, t_shift=calc_kl(par1.data[i],par2.data[k])
            if kl < min_kl:
                min_kl = kl
                j = k
                shift = t_shift

        if True: #(min_kl < kl_bound):
            common+=1

            # add logo pathts
            logo_path1 = os.path.join(
                logo_dir1,
                'logo{}.png'.format(i+1))
            logo_path2 = os.path.join(
                logo_dir2,
                'logo{}.png'.format(j+1))

            # add meta info
            meta_info1 = 'MOTIF {}: {} sequences'.format(i+1, par1.data[i]['nsites'])
            meta_info2 = 'MOTIF {}: {} sequences'.format(j+1, par2.data[j]['nsites'])

            logos.append((logo_path1, logo_path2))
            meta.append((meta_info1, meta_info2))

            print "{0:3}:{1:3}  {2:3}   {3:3}".format(i+1,j+1,min_kl,shift)
 
    print "\n"
    print "common number of motifs: {0}".format(common)
    
    return logos, meta

def generate_html(logos, meta):
    file_name = 'compare.html'
    global_meta_info = "This is an experiment on 11,642 sequences"
    meta_info1 = "Motifs returned by MEME"
    meta_info2 = "Motifs returned by ASC_MEME"
    
    doc, tag, text = Doc().tagtext()

    with tag('html'):
        with tag('body'):
            with tag('table', width = '500',  border="1", cellpadding="5"):
                
                # header text
                with tag('tr', align="center", valign="center"):
                    text(global_meta_info)
                with tag('tr'):
                    with tag('td', align="center", valign="center"):
                        text(meta_info1)

                    with tag('td', align="center", valign="center"):
                        text(meta_info2)
                
                for i in xrange(len(logos)):
                    logo_path1 = logos[i][0]
                    logo_path2 = logos[i][1]
                    meta_info1 = meta[i][0]
                    meta_info2 = meta[i][1]

                    with tag('tr'):
                        with tag('td', align="center", valign="center"):
                            doc.stag('img', src=logo_path1)
                            doc.stag('br')
                            text(meta_info1)

                        with tag('td', align="center", valign="center"):
                            doc.stag('img', src=logo_path2)
                            doc.stag('br')
                            text(meta_info2)
    
    with open(file_name, 'w') as f:
        f.write(doc.getvalue())
        
if __name__ == '__main__':
    '''
    [1] pwm file 1
    [2] pwm file 2
    '''
    out_dir1=sys.argv[1]
    out_dir2=sys.argv[2]
    
    logo_dir1 = os.path.join(out_dir1, 'logos')
    logo_dir2 = os.path.join(out_dir2, 'logos')
    
    pwm_file1 = os.path.join(out_dir1, 'pwms.txt')
    pwm_file2 = os.path.join(out_dir2, 'pwms.txt')
    
    par1=parse_motif.Parser([pwm_file1])
    par2=parse_motif.Parser([pwm_file2])
    
    logos, meta = compare(par1, logo_dir1, par2, logo_dir2)
    
    generate_html(logos, meta)

