import os
import sys
#import glob

from Queue import *
from threading import *

num_worker_threads = int(sys.argv[3])  # number of available workers (cores)
q = Queue()

def worker():
   while True:
       task = q.get()
       os.system(task)
       q.task_done() # important signal for q.join() to work




#guide_input_file_path = "../../../Guide_files_with_all_info/Gene_guides_files/"
guide_input_file_path = sys.argv[1]
script_to_be_executed = sys.argv[2]

all_input_guide_files = os.listdir(guide_input_file_path)


tasks = []

for files in all_input_guide_files:
   
   guide_input_file = os.path.join(guide_input_file_path, files)
   tasks.append("python %s  %s" %  (script_to_be_executed, guide_input_file))



for i in range(num_worker_threads):
    t = Thread(target=worker)
    t.daemon = True
    t.start()

for task in tasks:
   q.put(task)

q.join()
