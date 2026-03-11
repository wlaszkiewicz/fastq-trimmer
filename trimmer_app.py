import sys
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QScrollArea
from sections.plot_section import PlotsSection
from sections.trim_section    import TrimSection
from sections.derep_section   import DerepSection
from sections.cluster_section import ClusterSection
from sections.stats_section   import StatsSection
import os

class TrimmerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FASTQ Trimmer")
        self.setMinimumWidth(600)
        self.setMinimumHeight(600)
        self._build_ui()

    def _build_ui(self):
        self.trim_sec    = TrimSection()
        self.derep_sec   = DerepSection()
        self.cluster_sec = ClusterSection()
        self.plot_sec    = PlotsSection()
        self.stats_sec   = StatsSection()

        # wire autofill chain
        self.trim_sec.trim_finished.connect(self._on_trim_done)
        self.derep_sec.derep_finished.connect(self._on_derep_done)
        self.cluster_sec.cluster_finished.connect(self._on_cluster_done)


        # layout inside scroll area so it never gets squished
        inner = QWidget()
        inner_layout = QVBoxLayout()
        inner_layout.setSpacing(12)
        inner_layout.setContentsMargins(16, 16, 16, 16)
        inner_layout.addWidget(self.trim_sec)
        inner_layout.addWidget(self.derep_sec)
        inner_layout.addWidget(self.cluster_sec)
        inner_layout.addWidget(self.plot_sec)
        inner_layout.addWidget(self.stats_sec)
   
        inner_layout.addStretch()
        inner.setLayout(inner_layout)

        scroll = QScrollArea()
        scroll.setWidget(inner)
        scroll.setWidgetResizable(True)

        root = QVBoxLayout()
        root.setContentsMargins(0, 0, 0, 0)
        root.addWidget(scroll)
        self.setLayout(root)

    def _on_trim_done(self, trimmed_path, output_dir):
        self.derep_sec.autofill(trimmed_path, output_dir)
        self.stats_sec.autofill(trimmed_path, output_dir)

    def _on_derep_done(self, derep_path):
        output_dir = os.path.dirname(derep_path)
        self.cluster_sec.autofill(derep_path, output_dir)
        self.plot_sec.autofill(derep_path, output_dir)
   

    def _on_cluster_done(self, clustered_path, cluster_method):
        output_dir = os.path.dirname(clustered_path)
        method = "_size" if cluster_method == "cluster_size" else "_fast"
        self.plot_sec.autofill(f"clustered_sorted{method}.fasta", output_dir)
   

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = TrimmerApp()
    window.show()
    sys.exit(app.exec_())
