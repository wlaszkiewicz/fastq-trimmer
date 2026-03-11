import os
import re
from PyQt5.QtWidgets import (
    QGroupBox, QVBoxLayout, QLabel, QLineEdit, QPushButton,
    QCheckBox, QProgressBar, QHBoxLayout, QFormLayout, QMessageBox, QFileDialog
)
from PyQt5.QtCore import pyqtSignal
from trimmer.trimmer_worker import TrimWorker
from constants import ref1, ref2


class TrimSection(QGroupBox):
    # emits (trimmed_path, output_dir) when done
    trim_finished = pyqtSignal(str, str)

    def __init__(self, parent=None):
        super().__init__("1. Trim && Reverse", parent)
        self._build()

    def _build(self):
        layout = QVBoxLayout()
        layout.setSpacing(6)

        # --- primer sequences ---
        ref_group = QGroupBox("Primer sequences")
        ref_layout = QFormLayout()

        self.ref1_edit = QLineEdit(ref1)
        self.ref2_edit = QLineEdit(ref2)
        self.ref1_edit.setPlaceholderText("ref1 (forward start)")
        self.ref2_edit.setPlaceholderText("ref2 (forward end)")
        self.ref1_edit.setMaxLength(50)
        self.ref2_edit.setMaxLength(50)
        self.ref1_edit.setToolTip(
            "Forward start primer.\n"
            "Everything before this sequence will be cut off.\n"
            "The trimmed read will start here."
        )
        self.ref2_edit.setToolTip(
            "Forward end primer.\n"
            "Everything after this sequence will be cut off.\n"
            "The trimmed read will end here (inclusive)."
        )
        self.ref1_edit.textChanged.connect(
            lambda: self.ref1_edit.setText(self.ref1_edit.text().upper())
        )
        self.ref2_edit.textChanged.connect(
            lambda: self.ref2_edit.setText(self.ref2_edit.text().upper())
        )

        ref_layout.addRow("Forward start:", self.ref1_edit)
        ref_layout.addRow("Forward end:",   self.ref2_edit)
        ref_group.setLayout(ref_layout)
        layout.addWidget(ref_group)

        # --- barcode / input file ---
        self.barcode_label = QLabel("Barcode folder:")
        layout.addWidget(self.barcode_label)
        self.barcode_edit, row = self._path_row(
            "Select barcode folder...", self._pick_barcode, folder=True
        )
        layout.addLayout(row)

        # --- output folder ---
        layout.addWidget(QLabel("Output folder:"))
        self.output_edit, row = self._path_row(
            "Select output folder...", self._pick_output, folder=True
        )
        layout.addLayout(row)

        # --- merge checkbox ---
        self.merge_check = QCheckBox("Merge all .fastq.gz files before trimming")
        self.merge_check.setChecked(True)
        self.merge_check.stateChanged.connect(self._on_merge_toggle)
        layout.addWidget(self.merge_check)

        # --- status + progress ---
        self.status_lbl = QLabel("Ready.")
        layout.addWidget(self.status_lbl)
        self.progress_bar = QProgressBar()
        self.progress_bar.setValue(0)
        layout.addWidget(self.progress_bar)

        # --- run button ---
        self.run_btn = QPushButton("Trim && Reverse")
        self.run_btn.clicked.connect(self._run)
        layout.addWidget(self.run_btn)

        self.setLayout(layout)

    # helpers
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

    def _pick_barcode(self):
        if self.merge_check.isChecked():
            d = QFileDialog.getExistingDirectory(self, "Select barcode folder")
            if d:
                self.barcode_edit.setText(d)
        else:
            f, _ = QFileDialog.getOpenFileName(
                self, "Select fastq.gz file", "", "FASTQ GZ (*.fastq.gz *.gz)"
            )
            if f:
                self.barcode_edit.setText(f)

    def _pick_output(self):
        d = QFileDialog.getExistingDirectory(self, "Select output folder")
        if d:
            self.output_edit.setText(d)

    def _on_merge_toggle(self):
        if self.merge_check.isChecked():
            self.barcode_label.setText("Barcode folder:")
            self.barcode_edit.setPlaceholderText("Select barcode folder...")
        else:
            self.barcode_label.setText("Input file:")
            self.barcode_edit.setPlaceholderText("Select .fastq.gz file...")
        self.barcode_edit.clear()

    def _validate_ref(self, seq):
        return bool(re.fullmatch(r'[ACGTUacgtu]+', seq))

    def get_output_dir(self):
        return self.output_edit.text().strip()

    # run
    def _run(self):
        barcode_dir = self.barcode_edit.text().strip()
        output_dir  = self.output_edit.text().strip()
        ref1_val    = self.ref1_edit.text().strip().upper()
        ref2_val    = self.ref2_edit.text().strip().upper()

        if not ref1_val or not ref2_val:
            QMessageBox.warning(self, "Missing references", "Please enter both reference sequences.")
            return
        if not self._validate_ref(ref1_val) or not self._validate_ref(ref2_val):
            QMessageBox.warning(self, "Invalid references", "References can only contain A, C, G, T.")
            return
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
        self.progress_bar.setValue(0)

        self._worker = TrimWorker(
            barcode_dir, output_dir,
            self.merge_check.isChecked(),
            ref1_val, ref2_val
        )
        self._worker.progress.connect(self.progress_bar.setValue)
        self._worker.status.connect(self.status_lbl.setText)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_finished(self, complete, incomplete):
        self.run_btn.setEnabled(True)
        self.status_lbl.setText(f"Done! Complete: {complete}  Incomplete: {incomplete}")
        self.progress_bar.setValue(100)
        QMessageBox.information(
            self, "Trim done",
            f"Trimming complete!\n\nComplete:   {complete}\nIncomplete: {incomplete}"
        )
        trimmed_path = os.path.join(self.output_edit.text().strip(), "trimmed.fastq.gz")
        self.trim_finished.emit(trimmed_path, self.output_edit.text().strip())

    def _on_error(self, msg):
        self.run_btn.setEnabled(True)
        self.status_lbl.setText("Error.")
        QMessageBox.critical(self, "Trim error", msg)
