# Created by zty on 2018/10/24
import os
import glob
import argparse
import re
import shutil

process_dir = open("processed_dir","a+")

def parse_args():
    parser = argparse.ArgumentParser(description="保留对应的文件")
    parser.add_argument("-d",help="总目录")
    parser.add_argument("-ld",help="log文件夹目录")
    parser.add_argument("-k", metavar='N', type=str, nargs='+',help="保留文件")
    args = parser.parse_args()
    return args.d, args.ld, args.k

def clean_save(work_dir,log_dir,keep_files):
    os.chdir(work_dir)
    finished_samples = []
    for _file in glob.glob("%s/*success" %(log_dir)):
        sample_name = re.search(r".*\.(.*-[NT])\..*",_file).group(1)
        finished_samples.append(sample_name)
    for sample in finished_samples:
        os.chdir("Sample_%s" %(sample))
        process_dir.write(os.getcwd()+"\n")
        for _file in glob.glob("*"):
            if delete(_file,keep_files):
                pass
            else:
                try:
                    os.remove(_file)
                except Exception as e:
                    shutil.rmtree(_file, ignore_errors=True)
        os.chdir(work_dir)

def delete(_file,keep_files):
    for keep in keep_files:
        if _file.endswith(keep):
            return 1
        else:
            pass
    return

if __name__ == "__main__":
    work_dir,log_dir,keep_files = parse_args()
    clean_save(work_dir,log_dir,keep_files)
