import time
import psutil
from subprocess import PIPE
from concurrent.futures import ThreadPoolExecutor as pool
import threading
from threading import Timer


class Command():
    def __init__(self, cmd, name, timeout=3600*24*7):
        self.proc = None
        self.cmd = cmd
        self.name = name
        self.stderr = None
        self.stdout = None
        self.returncode = None
        self.timeout = timeout
        self.start_time = time.time()
        self.max_mem = 0
        self.max_cpu = 0

    def run(self):
        self.start_time = time.time()
        self.proc = psutil.Popen(self.cmd, shell=True, stderr=PIPE, stdout=PIPE)

        thread = threading.Thread(target=self.monitor_cmd, daemon=True)
        thread.start()

        timer = Timer(self.timeout, self.proc.kill)
        try:
            timer.start()
            self.stdout, self.stderr = self.proc.communicate()
        finally:
            timer.cancel()
        self.returncode = self.proc.returncode
        self.write_log()
        return self.returncode, self.max_mem, self.max_cpu

    def monitor_cmd(self):
        max_cpu = 0
        max_mem = 0
        while True:
            if self.proc is None:
                continue
            if self.proc.is_running():
                cpu_num = self.proc.cpu_num()
                if cpu_num > max_cpu:
                    max_cpu = cpu_num
                memory = self.proc.memory_info().rss / 1024 / 1024
                if memory > max_mem:
                    max_cpu = memory
            else:
                self.max_cpu = max_cpu
                self.max_mem = max_mem
                break
            time.sleep(1)
            if time.time() - self.start_time > self.timeout:
                self.max_cpu = max_cpu
                self.max_mem = max_mem
                break

    def write_log(self):
        if self.stderr:
            with open(self.name+'.'+str(self.proc.pid)+'.stderr', 'wb') as f:
                f.write(self.stderr)
        if self.stdout:
            with open(self.name+'.'+str(self.proc.pid)+'.stdout', 'wb') as f:
                f.write(self.stdout)


class Resource():
    def available_mem(self):
        pass

    def available_cpu(self):
        pass

    def is_enough(self, cpu=0, mem=0):
        if cpu >= self.available_cpu() and mem >= self.available_mem():
            return True
        else:
            return False

    def final_enough(self, cpu=0, mem=0, timeout=3600):
        start_time = time.time()
        while True:
            if self.is_enough(cpu, mem):
                return True
            if time.time() - start_time >= timeout:
                return None
            time.sleep(10)






cmd = Command('sleep 5 && echo hello-world', 'sleep')
print(cmd.run())
