import sys
import os
from PyQt5.QtWidgets import (
    QApplication, QWidget, QLabel, QPushButton, QLineEdit,
    QFileDialog, QCheckBox, QProgressBar, QVBoxLayout,
    QHBoxLayout, QGroupBox, QMessageBox, QDialog, QVBoxLayout, QDoubleSpinBox
)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from trimmer_worker import TrimWorker
from vsearch_worker import DerepWorker, StatsWorker
from plot_logic import plot_length_distribution


class TrimmerApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("FASTQ Trimmer")
        self.setMinimumWidth(540)
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout()
        root.setSpacing(10)
        root.setContentsMargins(16, 16, 16, 16)

        root.addWidget(self._build_trim_section())
        root.addWidget(self._build_vsearch_section())

        self.setLayout(root)

    def _build_trim_section(self):
        box = QGroupBox("1. Trim && Reverse")
        layout = QVBoxLayout()
        layout.setSpacing(6)

        # barcode folder
        self.barcode_label = QLabel("Barcode folder:")
        layout.addWidget(self.barcode_label)
        self.barcode_edit, row = self._path_row("Select barcode folder...", self._pick_barcode, folder=True)
        layout.addLayout(row)

        # output folder
        layout.addWidget(QLabel("Output folder:"))
        self.output_edit, row = self._path_row("Select output folder...", self._pick_output, folder=True)
        self.output_edit.textChanged.connect(self._autofill_vsearch_input)
        layout.addLayout(row)

        # merge checkbox
        self.merge_check = QCheckBox("Merge all .fastq.gz files before trimming")
        self.merge_check.setChecked(True)
        layout.addWidget(self.merge_check)
        self.merge_check.stateChanged.connect(self._on_merge_toggle)


        # status + progress
        self.trim_status = QLabel("Ready.")
        layout.addWidget(self.trim_status)
        self.trim_progress = QProgressBar()
        self.trim_progress.setValue(0)
        layout.addWidget(self.trim_progress)

        # run button
        self.run_btn = QPushButton("Trim && Reverse")
        self.run_btn.clicked.connect(self._run_trim)
        layout.addWidget(self.run_btn)

        box.setLayout(layout)
        return box

    def _build_vsearch_section(self):
        box = QGroupBox("2. vsearch")
        layout = QVBoxLayout()
        layout.setSpacing(6)

        layout.addWidget(QLabel("Input file (trimmed.fastq.gz):"))
        self.vsearch_input_edit, row = self._path_row("Auto-filled after trim, or browse...", self._pick_vsearch_input, folder=False)
        layout.addLayout(row)

        # vsearch output folder
        layout.addWidget(QLabel("Output folder:"))
        self.vsearch_output_edit, row = self._path_row("Select output folder...", self._pick_vsearch_output, folder=True)
        layout.addLayout(row)

        # status
        self.vsearch_status = QLabel("Ready.")
        layout.addWidget(self.vsearch_status)

        # buttons row
        btn_row = QHBoxLayout()
        self.derep_btn = QPushButton("Run Derep")
        self.derep_btn.setToolTip(
            "Merge strictly identical sequences contained in a FASTQ file.\n"
            "Identical sequences are defined as having the same length and the same\n"
            "string of nucleotides (case insensitive, T and U are considered the same)."
        )
        self.derep_btn.clicked.connect(self._run_derep)
        self.stats_btn = QPushButton("Run Stats")
        self.stats_btn.setToolTip(
            "Analyze fastq sequences, providing 5 detailed statistics tables.\n"
            "Including read length distribution, quality score distribution, and more.\n"
            "Learn more at: https://torognes.github.io/vsearch/commands/vsearch-fastq_stats.1.html"
        )
        self.stats_btn.clicked.connect(self._run_stats)
        btn_row.addWidget(self.derep_btn)
        btn_row.addWidget(self.stats_btn)
        layout.addLayout(btn_row)

        plot_row = QHBoxLayout()
        self.plot_btn = QPushButton("Plot Length Distribution")
        self.plot_btn.clicked.connect(self._run_plot)
        plot_row.addWidget(self.plot_btn, stretch=1)

        # right half: threshold centered inside it
        right_half = QHBoxLayout()
        right_half.addStretch()
        right_half.addWidget(QLabel("Outlier threshold (%):"))
        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(0.0, 10.0)
        self.threshold_spin.setSingleStep(0.05)
        self.threshold_spin.setValue(0.1)
        self.threshold_spin.setDecimals(2)
        self.threshold_spin.setFixedWidth(70)
        self.threshold_spin.setToolTip(
            "Hide lengths whose read count is below this % of total reads.\n"
            "Higher = more aggressive cropping. 0.0 = show everything."
        )
        right_half.addWidget(self.threshold_spin)
        right_half.addStretch()

        right_widget = QWidget()
        right_widget.setLayout(right_half)
        plot_row.addWidget(right_widget, stretch=1)
        layout.addLayout(plot_row)

        box.setLayout(layout)
        return box

    #path row helper 
    def _path_row(self, placeholder, callback, folder=True):
        edit = QLineEdit()
        edit.setPlaceholderText(placeholder)
        edit.setReadOnly(True)
        btn = QPushButton("Browse")
        btn.setFixedWidth(70)
        btn.clicked.connect(callback)
        row = QHBoxLayout()
        row.addWidget(edit)
        row.addWidget(btn)
        return edit, row

    # folder / file pickers
    def _pick_barcode(self):
        if self.merge_check.isChecked():
            d = QFileDialog.getExistingDirectory(self, "Select barcode folder")
            if d:
                self.barcode_edit.setText(d)
        else:
            f, _ = QFileDialog.getOpenFileName(self, "Select fastq.gz file", "", "FASTQ GZ (*.fastq.gz *.gz)")
            if f:
                self.barcode_edit.setText(f)

    def _pick_output(self):
        d = QFileDialog.getExistingDirectory(self, "Select output folder")
        if d:
            self.output_edit.setText(d)

    def _pick_vsearch_input(self):
        f, _ = QFileDialog.getOpenFileName(self, "Select input fastq.gz", "", "FASTQ GZ (*.fastq.gz *.gz)")
        if f:
            self.vsearch_input_edit.setText(f)

    def _pick_vsearch_output(self):
        d = QFileDialog.getExistingDirectory(self, "Select vsearch output folder")
        if d:
            self.vsearch_output_edit.setText(d)

    def _autofill_vsearch_input(self, output_dir):
        """When output folder is set, pre-fill vsearch input with trimmed.fastq.gz path."""
        trimmed_path = os.path.join(output_dir, "trimmed.fastq.gz")
        self.vsearch_input_edit.setText(trimmed_path)
        self.vsearch_output_edit.setText(output_dir)
    
    def _on_merge_toggle(self):
        if self.merge_check.isChecked():
            self.barcode_label.setText("Barcode folder:")
            self.barcode_edit.setPlaceholderText("Select barcode folder...")
        else:
            self.barcode_label.setText("Input file:")
            self.barcode_edit.setPlaceholderText("Select .fastq.gz file...")
        self.barcode_edit.clear()

    def _run_plot(self):
        log_path   = os.path.join(self.vsearch_output_edit.text().strip(), "fastq_stats.log")
        output_dir = self.vsearch_output_edit.text().strip()

        if not os.path.isfile(log_path):
            QMessageBox.warning(self, "Missing log", f"Could not find fastq_stats.log in the output folder.\nRun Stats first!")
            return

        try:
            fig = plot_length_distribution(log_path, threshold_pct=round(self.threshold_spin.value(), 2))

            dialog = QDialog(self)
            dialog.setWindowTitle("Read Length Distribution")
            dialog.resize(900, 450)
            layout = QVBoxLayout()
            canvas = FigureCanvasQTAgg(fig)
            layout.addWidget(canvas)
            btn_row = QHBoxLayout()

            save_btn = QPushButton("Save as PNG")
            close_btn = QPushButton("Close")

            def save_plot():
                output_png = os.path.join(output_dir, "length_distribution.png")
                fig.savefig(output_png, dpi=150)
                QMessageBox.information(dialog, "Saved", f"Plot saved to:\n{output_png}")

            save_btn.clicked.connect(save_plot)
            close_btn.clicked.connect(dialog.accept)

            btn_row.addWidget(save_btn)
            btn_row.addStretch()
            btn_row.addWidget(close_btn)
            layout.addLayout(btn_row)

            dialog.setLayout(layout)
            dialog.exec_()

            plt.close(fig)

        except Exception as e:
            QMessageBox.critical(self, "Plot error", str(e))

    # trim
    def _run_trim(self):
        barcode_dir = self.barcode_edit.text().strip()
        output_dir  = self.output_edit.text().strip()

        if self.merge_check.isChecked():
            if not barcode_dir or not os.path.isdir(barcode_dir):
                QMessageBox.warning(self, "Missing input", "Please select a valid barcode folder.")
                return
        else:
            if not barcode_dir or not os.path.isfile(barcode_dir):
                QMessageBox.warning(self, "Missing input", "Please select a valid .fastq.gz file.")
                return

        if not output_dir or not os.path.isdir(output_dir):
            QMessageBox.warning(self, "Missing output", "Please select a valid output folder.")
            return

        self.run_btn.setEnabled(False)
        self.trim_progress.setValue(0)

        self.trim_worker = TrimWorker(barcode_dir, output_dir, self.merge_check.isChecked())
        self.trim_worker.progress.connect(self.trim_progress.setValue)
        self.trim_worker.status.connect(self.trim_status.setText)
        self.trim_worker.finished.connect(self._on_trim_finished)
        self.trim_worker.error.connect(self._on_trim_error)
        self.trim_worker.start()

    def _on_trim_finished(self, complete, incomplete):
        self.run_btn.setEnabled(True)
        self.trim_status.setText(f"Done! Complete: {complete}  Incomplete: {incomplete}")
        self.trim_progress.setValue(100)
        QMessageBox.information(
            self, "Trim done",
            f"Trimming complete!\n\nComplete:   {complete}\nIncomplete: {incomplete}\n\nFiles saved to output folder."
        )

    def _on_trim_error(self, msg):
        self.run_btn.setEnabled(True)
        self.trim_status.setText("Error.")
        QMessageBox.critical(self, "Trim error", msg)

    # derep 
    def _run_derep(self):
        input_path = self.vsearch_input_edit.text().strip()
        output_dir = self.vsearch_output_edit.text().strip()

        if not input_path or not os.path.isfile(input_path):
            QMessageBox.warning(self, "Missing input", "Please select a valid input .fastq.gz file.")
            return
        if not output_dir or not os.path.isdir(output_dir):
            QMessageBox.warning(self, "Missing output", "Please select a valid output folder.")
            return

        self.derep_btn.setEnabled(False)
        self.stats_btn.setEnabled(False)
        self.plot_btn.setEnabled(False)
        self.vsearch_status.setText("Starting derep...")

        self.derep_worker = DerepWorker(input_path, output_dir)
        self.derep_worker.status.connect(self.vsearch_status.setText)
        self.derep_worker.finished.connect(self._on_derep_finished)
        self.derep_worker.error.connect(self._on_vsearch_error)
        self.derep_worker.start()

    def _on_derep_finished(self, output_path):
        self.derep_btn.setEnabled(True)
        self.stats_btn.setEnabled(True)
        self.plot_btn.setEnabled(True)
        QMessageBox.information(self, "Derep done", f"Dereplication complete!\n\nOutput: {output_path}")

    # stats 
    def _run_stats(self):
        input_path = self.vsearch_input_edit.text().strip()
        output_dir = self.vsearch_output_edit.text().strip()

        if not input_path or not os.path.isfile(input_path):
            QMessageBox.warning(self, "Missing input", "Please select a valid input .fastq.gz file.")
            return
        if not output_dir or not os.path.isdir(output_dir):
            QMessageBox.warning(self, "Missing output", "Please select a valid output folder.")
            return

        self.derep_btn.setEnabled(False)
        self.stats_btn.setEnabled(False)
        self.plot_btn.setEnabled(False)
        self.vsearch_status.setText("Starting stats...")

        self.stats_worker = StatsWorker(input_path, output_dir)
        self.stats_worker.status.connect(self.vsearch_status.setText)
        self.stats_worker.finished.connect(self._on_stats_finished)
        self.stats_worker.error.connect(self._on_vsearch_error)
        self.stats_worker.start()

    def _on_stats_finished(self, log_path):
        self.derep_btn.setEnabled(True)
        self.stats_btn.setEnabled(True)
        self.plot_btn.setEnabled(True)
        QMessageBox.information(self, "Stats done", f"Stats complete!\n\nLog saved to: {log_path}")

    def _on_vsearch_error(self, msg):
        self.derep_btn.setEnabled(True)
        self.stats_btn.setEnabled(True)
        self.plot_btn.setEnabled(True)
        self.vsearch_status.setText("Error.")
        QMessageBox.critical(self, "vsearch error", msg)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = TrimmerApp()
    window.show()
    sys.exit(app.exec_())
