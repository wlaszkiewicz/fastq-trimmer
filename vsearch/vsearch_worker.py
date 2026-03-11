import os
from PyQt5.QtCore import QThread, pyqtSignal
from vsearch.vsearch_logic import run_derep, run_fastq_stats, run_clustering

class DerepWorker(QThread):
    status   = pyqtSignal(str)
    finished = pyqtSignal(str)   
    error    = pyqtSignal(str)

    def __init__(self, input_path, output_dir):
        super().__init__()
        self.input_path = input_path
        self.output_dir = output_dir

    def run(self):
        try:
            # run derep 
            self.status.emit("Running vsearch --fastx_uniques...")
            derep_path = os.path.join(self.output_dir, "derep.fasta")
            run_derep(self.input_path, derep_path)

            self.status.emit(f"Done. Output: {derep_path}")
            self.finished.emit(derep_path)

        except Exception as e:
            self.error.emit(str(e))


class StatsWorker(QThread):
    status   = pyqtSignal(str)
    finished = pyqtSignal(str)  
    error    = pyqtSignal(str)

    def __init__(self, input_path, output_dir):
        super().__init__()
        self.input_path = input_path
        self.output_dir = output_dir

    def run(self):
        try:
            self.status.emit("Running vsearch --fastq_stats...")
            log_path = os.path.join(self.output_dir, "fastq_stats.log")
            run_fastq_stats(self.input_path, log_path)
            self.status.emit(f"Done. Log saved: {log_path}")
            self.finished.emit(log_path)

        except Exception as e:
            self.error.emit(str(e))


class ClusterWorker(QThread):
    status   = pyqtSignal(str)
    finished = pyqtSignal(str)  
    error    = pyqtSignal(str)

    def __init__(self, input_path, output_dir, method='cluster_size', identity=0.97, minsize=1):
        super().__init__()
        self.input_path = input_path
        self.output_dir = output_dir
        self.type       = method
        self.identity   = identity
        self.minsize    = minsize

    def run(self):
        try:
            self.status.emit("Running vsearch clustering...")
            uc_path     = os.path.join(self.output_dir, "clusters.uc")
            sorted_path = run_clustering(
                self.input_path, self.output_dir,
                self.type, self.identity, self.minsize,
                status_callback=self.status.emit
            )
            self.status.emit(f"Done. Output: clusters: {uc_path}, sorted: {sorted_path}")
            self.finished.emit(sorted_path)

        except Exception as e:
            self.error.emit(str(e))


