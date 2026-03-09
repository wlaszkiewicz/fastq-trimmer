import os
from PyQt5.QtCore import QThread, pyqtSignal
from trimmer_logic import merge_fastq_gz, process_file


class TrimWorker(QThread):
    progress = pyqtSignal(int)       # 0-100
    status   = pyqtSignal(str)       # status message
    finished = pyqtSignal(int, int)  # complete, incomplete
    error    = pyqtSignal(str)

    def __init__(self, barcode_dir, output_dir, do_merge):
        super().__init__()
        self.barcode_dir = barcode_dir
        self.output_dir  = output_dir
        self.do_merge    = do_merge

    def run(self):
        try:
            # merge step 
            if self.do_merge:
                self.status.emit("Merging .fastq.gz files...")
                merged_path = os.path.join(self.output_dir, "merged.fastq.gz")
                merge_fastq_gz(self.barcode_dir, merged_path)
            else:
                merged_path = self.barcode_dir

            # trim step
            self.status.emit("Counting records...")
            output_path     = os.path.join(self.output_dir, "trimmed.fastq.gz")
            incomplete_path = os.path.join(self.output_dir, "incomplete.fastq.gz")

            self.status.emit("Trimming...")
            complete, incomplete = process_file(
                merged_path, output_path, incomplete_path,
                progress_callback=self.progress.emit
            )

            self.finished.emit(complete, incomplete)

        except Exception as e:
            self.error.emit(str(e))
