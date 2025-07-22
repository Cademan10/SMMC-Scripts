import os

os.chdir("..")

home = os.getcwd()
dirs = os.listdir(home)

dirs = [dir for dir in dirs if dir[0]=="B"]
run = input("Run all or specific beta?: ")
if run.lower()=="all":

    for dir in dirs:
        os.chdir(dir+"/dB64")
        os.system("sbatch run_calc.sh")

        os.chdir("../dB32")
        os.system("sbatch run_calc.sh")

        os.chdir(home)

else:
    dir = "Beta"+run
    os.chdir(dir+"/dB64")
    os.system("sbatch run_calc.sh")

    os.chdir("../dB32")
    os.system("sbatch run_calc.sh")

    os.chdir(home)

        
