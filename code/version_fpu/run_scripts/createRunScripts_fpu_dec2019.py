##################################################
# create run scripts for cell level optimization
##################################################

import os, stat
import platform


theOS= platform.system()

if(theOS== 'Windows'):
    basePath = "c:/sven/Nutzer/projekte_neu/ess_tradeoffs/LPJ_pareto/"
else:
    basePath = "/home/s/slautenbach/projects/optim/LPJ_pareto/"

def writeSubmitScript (fn, content):
    submitfile = open(fn, 'wb')
    submitfile.write(content)
    submitfile.close()

mutation_rate= 0.1
tournament_size= 2
max_generations= 300
pop_extra= 500
corFac= 0.98

codeFolder= basePath + "code/version_fpu_updated_dec_2019"

for runID in [21,22]:
    for aRcp in ["26", "60"]:
        for aTimePeriod in["2042", "2090"]:
            if aTimePeriod == "2042":
                aTimePeriodLong= "2033-2042"
            else:
                aTimePeriodLong = "2090-2099"
            content= "#!/bin/bash\n"
            content= content + "source ~/venv/optim2updt/bin/activate\n"
            content= content + "cd " + codeFolder + "\n\n"
            content= content + "python2.7 lpf_pareto_fpu_forage.py -r " + str(runID) + " --mutation_rate " + str(mutation_rate)
            content= content + " --tournament_size=" + str(tournament_size) + " --max_generations=" + str(max_generations)
            content= content + " --time_period=" + str(aTimePeriod) + " --rcp=" + str(aRcp) + " --pop_extra=" + str(pop_extra)
            content= content + " --corFac=" + str(corFac)  + "\n"

            fn= codeFolder + "/run_scripts/run_" + str(runID) + "_" + str(aRcp) + "_" + str(aTimePeriodLong) + ".sh"
            writeSubmitScript(fn, content)
            os.chmod(fn, stat.S_IRWXU)
