#
#11/20/2014 2:05:54 PM
#Author: Honglei Liu <liuhonglei@gmail.com>
#
#Description: read and parse results of MEME
#

import sys
import re

class Parser:
    "'extract the pwms from the given file'"
    #beta=0.0001
    def __init__(self,files):
        pattern=re.compile(r"^letter-probability matrix:\s+(.*)$")
        self.data=[]
        self.count=0
        for fname in files:
            with open(fname) as f:
                while True:
                    line=f.readline()
                    if len(line)==0:
                        break
                    if pattern.match(line) is not None:
                        self.data.append({})
                        sp=line.split()
                        self.data[self.count]['alength']=int(sp[3])
                        self.data[self.count]['w']=int(sp[5])
                        self.data[self.count]['nsites']=int(sp[7])
                        self.data[self.count]['evalue']=float(sp[9])

                        self.data[self.count]['pwm_p']=[]
                        self.data[self.count]['pwm_f']=[]                        
                        for i in range(self.data[self.count]['w']):
                            line=f.readline()
                            temp_p=[]
                            temp_f=[]
                            for lp in line.split():
                                t_p=float(lp)
                                #tp=(self.data[self.count]['nsites']*t+self.beta)/(self.data[self.count]['nsites']+self.beta*self.data[self.count]['alength'])
                                temp_p.append(t_p)
                                temp_f.append(self.data[self.count]['nsites']*t_p+5)  #at least have freq. 5
                            self.data[self.count]['pwm_p'].append(temp_p[:])
                            self.data[self.count]['pwm_f'].append(temp_f[:])
                        self.count+=1

if __name__ == '__main__':

    par=Parser(sys.argv[1:])
    print par.count
    print par.data[1]['evalue']
