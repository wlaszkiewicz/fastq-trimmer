import os
from PyQt5.QtCore import QThread, pyqtSignal
from vsearch_logic import run_derep, run_fastq_stats


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
