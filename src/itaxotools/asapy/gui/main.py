#-----------------------------------------------------------------------------
# ASAPy - Assemble Species by Automatic Partitioning with ASAP
# Copyright (C) 2021  Patmanidis Stefanos
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------


"""
GUI for ASAP ...
"""

from PySide6 import QtCore
from PySide6 import QtWidgets
from PySide6 import QtGui
from PySide6 import QtSvgWidgets
from PySide6 import QtStateMachine

import os
import sys
import logging
import tempfile
import time
import shutil
import pathlib
import re

from itaxotools.common.param.view import View as ParamView
from itaxotools.common.param.model import Model as ParamModel
from itaxotools.common import threading
from itaxotools.common import machine
from itaxotools.common import widgets
from itaxotools import common

from .. import core

from . import utility
from . import widgets
from . import resources


def work_run(analysis):
    analysis.run()
    return analysis.results


class ResultItem(QtWidgets.QListWidgetItem):
    """
    Model for an ASAP analysis result file.
    Holds icon, label, tooltip, category and file location.
    """
    Type = QtWidgets.QListWidgetItem.UserType + 1
    Icons = { None: QtGui.QIcon() }

    def __init__(self, file):
        """Overloaded with new type"""
        super().__init__(type=self.Type)
        self.file = file
        path = pathlib.Path(file)
        suffix = path.suffix
        if not (suffix in self.Icons.keys()):
            suffix = None
        self.setIcon(self.Icons[suffix])
        self.setText(path.name)

class ResultView(QtWidgets.QListWidget):
    """
    Show each result file with icon and label.
    Remember file type and path
    """
    def open(self, folder):
        """Refresh contents"""
        self.clear()
        path = pathlib.Path(folder)

        # result files
        for file in sorted(list(path.glob('*.tab')), reverse=True):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # graph files
        for file in sorted(list(path.glob('*.svg'))):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # spart files
        for file in sorted(list(path.glob('*.spart'))):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # spart XML files
        for file in sorted(list(path.glob('*.spart.xml'))):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # partition files (lists)
        for file in sorted(list(path.glob('*.txt'))):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # partition files (csv)
        for file in sorted(list(path.glob('*.csv'))):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # tree files
        for file in sorted(list(path.glob('*.tree'))):
            item = ResultItem(str(path / file))
            self.addItem(item)

        # log files
        for file in list(path.glob('*.log')):
            item = ResultItem(str(path / file))
            self.addItem(item)


class Main(common.widgets.ToolDialog):
    """Main window, handles everything"""

    actionSignal = QtCore.Signal(str, list, dict)

    def __init__(self, parent=None, init=None):
        super(Main, self).__init__(parent)

        self.logger = logging.getLogger('iTaxotools.ASAPy')
        self.logger.setLevel(logging.DEBUG)
        handler = logging.StreamHandler(sys.stdout)
        self.logger.addHandler(handler)

        self.title = 'ASAPy'
        self.analysis = core.PartitionAnalysis(None)
        self._temp = None
        self.temp = None

        self.setWindowTitle(self.title)
        self.setWindowIcon(QtGui.QIcon(resources.get('asap-icon-transparent.ico')))
        self.resize(1024,480)

        self.machine = None
        self.skin()
        self.draw()
        self.act()
        self.cog()

        self.setParamModel(self.analysis)

        if init is not None:
            self.machine.started.connect(lambda: init(self))

    def __getstate__(self):
        return (self.analysis,)

    def __setstate__(self, state):
        (self.analysis,) = state

    def setParamModel(self, analysis):
        self.param_model = ParamModel(analysis.params, self)
        self.param_model.dataChanged.connect(lambda x: self.postAction('UPDATE'))
        self.paramWidget.setModel(self.param_model)

    def postAction(self, action, *args, **kwargs):
        self.actionSignal.emit(action, args, kwargs)

    def taggedTransition(self, action):
        return machine.TaggedTransition(self.actionSignal, action)

    def filterReject(self):
        """If running, verify cancel"""
        if self.state['running'] in list(self.machine.configuration()):
            self.handleStop()
            return False
        else:
            return True

    def skin(self):
        """Configure widget appearance"""
        color = {
            'white':  '#ffffff',
            'light':  '#eff1ee',
            'beige':  '#e1e0de',
            'gray':   '#abaaa8',
            'iron':   '#8b8d8a',
            'black':  '#454241',
            'red':    '#ee4e5f',
            'pink':   '#eb9597',
            'orange': '#eb6a4a',
            'brown':  '#655c5d',
            'green':  '#00ff00',
            }
        # using green for debugging
        palette = QtGui.QGuiApplication.palette()
        scheme = {
            QtGui.QPalette.Active: {
                QtGui.QPalette.Window: 'light',
                QtGui.QPalette.WindowText: 'black',
                QtGui.QPalette.Base: 'white',
                QtGui.QPalette.AlternateBase: 'light',
                QtGui.QPalette.PlaceholderText: 'brown',
                QtGui.QPalette.Text: 'black',
                QtGui.QPalette.Button: 'light',
                QtGui.QPalette.ButtonText: 'black',
                QtGui.QPalette.Light: 'white',
                QtGui.QPalette.Midlight: 'beige',
                QtGui.QPalette.Mid: 'gray',
                QtGui.QPalette.Dark: 'iron',
                QtGui.QPalette.Shadow: 'brown',
                QtGui.QPalette.Highlight: 'red',
                QtGui.QPalette.HighlightedText: 'white',
                # These work on linux only?
                QtGui.QPalette.ToolTipBase: 'beige',
                QtGui.QPalette.ToolTipText: 'brown',
                # These seem bugged anyway
                QtGui.QPalette.BrightText: 'green',
                QtGui.QPalette.Link: 'green',
                QtGui.QPalette.LinkVisited: 'green',
                },
            QtGui.QPalette.Disabled: {
                QtGui.QPalette.Window: 'light',
                QtGui.QPalette.WindowText: 'iron',
                QtGui.QPalette.Base: 'white',
                QtGui.QPalette.AlternateBase: 'light',
                QtGui.QPalette.PlaceholderText: 'green',
                QtGui.QPalette.Text: 'iron',
                QtGui.QPalette.Button: 'light',
                QtGui.QPalette.ButtonText: 'gray',
                QtGui.QPalette.Light: 'white',
                QtGui.QPalette.Midlight: 'beige',
                QtGui.QPalette.Mid: 'gray',
                QtGui.QPalette.Dark: 'iron',
                QtGui.QPalette.Shadow: 'brown',
                QtGui.QPalette.Highlight: 'pink',
                QtGui.QPalette.HighlightedText: 'white',
                # These seem bugged anyway
                QtGui.QPalette.BrightText: 'green',
                QtGui.QPalette.ToolTipBase: 'green',
                QtGui.QPalette.ToolTipText: 'green',
                QtGui.QPalette.Link: 'green',
                QtGui.QPalette.LinkVisited: 'green',
                },
            }
        scheme[QtGui.QPalette.Inactive] = scheme[QtGui.QPalette.Active]
        for group in scheme:
            for role in scheme[group]:
                palette.setColor(group, role,
                    QtGui.QColor(color[scheme[group][role]]))
        QtGui.QGuiApplication.setPalette(palette)

        self.colormap = {
            widgets.VectorIcon.Normal: {
                '#000': color['black'],
                '#f00': color['red'],
                },
            widgets.VectorIcon.Disabled: {
                '#000': color['gray'],
                '#f00': color['orange'],
                },
            }
        self.colormap_icon =  {
            '#000': color['black'],
            '#ff0000': color['red'],
            '#ffa500': color['pink'],
            }
        self.colormap_icon_light =  {
            '#000': color['iron'],
            '#ff0000': color['red'],
            '#ffa500': color['pink'],
            }
        self.colormap_graph =  {
            'groups': {
                'black':   color['black'],
                },
            'histogram': {
                'black':   color['black'],
                '#ffb83b': color['beige'],
                },
            'ranks': {
                'black':   color['black'],
                '#00a5eb': color['red'],
                },
            'species': {
                'black':   color['black'],
                },
            }

        ResultItem.Icons['.tab'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-text.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.txt'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-text.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.csv'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-text.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.svg'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-graph.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.log'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-log.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.spart'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-spart.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.xml'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-spart.svg'),
                colormap=self.colormap_icon))
        ResultItem.Icons['.tree'] = \
            QtGui.QIcon(widgets.VectorPixmap(resources.get('file-tree.svg'),
                colormap=self.colormap_icon))

    def draw(self):
        """Draw all widgets"""
        self.header = widgets.Header()
        self.header.logoTool = widgets.VectorPixmap(resources.get('logo-asap.svg'),
            colormap=self.colormap_icon)
        self.header.logoProject = QtGui.QPixmap(resources.get('itaxotools-logo-new.png'))
        self.header.description = (
            'Primary species delimitation' + '\n'
            'using automatic partitioning'
            )
        self.header.citation = (
            'ASAP by S. Brouillet, G. Achaz, N. Puillandre' + '\n'
            'BIONJ by O. Gascuel, GUI by S. Patmanidis'
        )

        self.line = widgets.Subheader()

        self.line.icon = QtWidgets.QLabel()
        self.line.icon.setPixmap(widgets.VectorPixmap(resources.get('arrow-right.svg'),
            colormap=self.colormap_icon_light))
        self.line.icon.setStyleSheet('border-style: none;')

        self.line.file = QtWidgets.QLineEdit()
        self.line.file.setPlaceholderText('Open a file to begin')
        self.line.file.setReadOnly(True)
        self.line.file.setStyleSheet("""
            QLineEdit {
                background-color: palette(Base);
                padding: 2px 4px 2px 4px;
                border-radius: 4px;
                border: 1px solid palette(Mid);
                }
            """)

        layout = QtWidgets.QHBoxLayout()
        layout.addSpacing(4)
        layout.addWidget(self.line.icon)
        layout.addSpacing(2)
        layout.addWidget(self.line.file)
        layout.addSpacing(14)
        layout.setContentsMargins(4, 4, 4, 4)
        self.line.setLayout(layout)

        self.pane = {}

        self.paramWidget = ParamView()
        self.paramWidget.setSizePolicy(
            QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.MinimumExpanding)
        self.paramWidget.setContentsMargins(0, 0, 0, 0)

        self.pane['param'] = widgets.Panel(self)
        self.pane['param'].title = 'Parameters'
        self.pane['param'].footer = 'Hover parameters for tips'
        self.pane['param'].body.addWidget(self.paramWidget)
        self.pane['param'].body.addStretch(1)
        self.pane['param'].body.setContentsMargins(0, 0, 0, 0)

        self.folder = ResultView()
        self.folder.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.folder.itemActivated.connect(self.handlePreview)
        self.folder.setStyleSheet("ResultView::item {padding: 2px;}")
        self.folder.setAlternatingRowColors(True)

        self.pane['list'] = widgets.Panel(self)
        self.pane['list'].title = 'Files'
        self.pane['list'].footer = 'Nothing to show'
        self.pane['list'].body.addWidget(self.folder)
        self.pane['list'].body.setContentsMargins(1, 1, 1, 1)
        # self.pane['list'].body.addStretch(1)

        self.preview = QtWidgets.QTextEdit()
        self.preview.setFont(
            QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.preview.setWordWrapMode(QtGui.QTextOption.WrapMode(0))
        self.preview.setReadOnly(True)

        self.graphSvg = QtSvgWidgets.QSvgWidget()

        self.graph = QtWidgets.QFrame()
        self.graph.setStyleSheet("""
            QFrame {
                background-color: palette(Base);
                border: 1px solid palette(Mid);
                }
            """)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.graphSvg)
        self.graph.setLayout(layout)

        self.stack = QtWidgets.QStackedLayout()
        self.stack.addWidget(self.preview)
        self.stack.addWidget(self.graph)

        self.pane['preview'] = widgets.Panel(self)
        self.pane['preview'].title = 'Preview'
        self.pane['preview'].footer = 'Nothing to show'
        self.pane['preview'].body.addLayout(self.stack)
        self.pane['preview'].body.setContentsMargins(1, 1, 1, 1)

        self.splitter = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        self.splitter.addWidget(self.pane['param'])
        self.splitter.addWidget(self.pane['list'])
        self.splitter.addWidget(self.pane['preview'])
        self.splitter.setStretchFactor(0,0)
        self.splitter.setStretchFactor(1,0)
        self.splitter.setStretchFactor(2,1)
        self.splitter.setCollapsible(0,False)
        self.splitter.setCollapsible(1,False)
        self.splitter.setCollapsible(2,False)
        self.splitter.setStyleSheet("QSplitter::handle { height: 12px; }")
        self.splitter.setContentsMargins(8, 4, 8, 4)
        self.splitter.setSizes([
            self.width()/4,
            self.width()/4,
            self.width()/2,
            ])

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.header)
        layout.addWidget(self.line)
        layout.addSpacing(8)
        layout.addWidget(self.splitter)

        layout.setSpacing(0)
        layout.setContentsMargins(0, 0, 0, 0)
        self.setLayout(layout)

        self.setContentsMargins(0, 0, 0, 0)

    def act(self):
        """Populate dialog actions"""
        self.action = {}

        self.action['open'] = QtGui.QAction('&Open', self)
        self.action['open'].setIcon(widgets.VectorIcon(resources.get('open.svg'), self.colormap))
        self.action['open'].setShortcut(QtGui.QKeySequence.Open)
        self.action['open'].setToolTip((
            'Open an aligned fasta file or a distance matrix\n'
            'Accepted formats: phylip, dnadist and MEGA'))
        self.action['open'].triggered.connect(self.handleOpen)

        self.action['save'] = QtGui.QAction('&Save', self)
        self.action['save'].setIcon(widgets.VectorIcon(resources.get('save.svg'), self.colormap))
        self.action['save'].setShortcut(QtGui.QKeySequence.Save)
        self.action['save'].setToolTip((
            'Save files with a prefix of your choice\n'
            'Change filter to choose what files are saved'))
        self.action['save'].triggered.connect(self.handleSave)

        self.action['run'] = QtGui.QAction('&Run', self)
        self.action['run'].setIcon(widgets.VectorIcon(resources.get('run.svg'), self.colormap))
        self.action['run'].setShortcut('Ctrl+R')
        self.action['run'].setToolTip('Run ASAP analysis')
        self.action['run'].triggered.connect(self.handleRun)

        self.action['stop'] = QtGui.QAction('Stop', self)
        self.action['stop'].setIcon(widgets.VectorIcon(resources.get('stop.svg'), self.colormap))
        # self.action['stop'].setShortcut(QtGui.QKeySequence.Cancel)
        self.action['stop'].setToolTip('Cancel analysis')
        self.action['stop'].triggered.connect(self.handleStop)

        for action in self.action.values():
            self.header.toolbar.addAction(action)

    def cog(self):
        """Define state machine"""

        self.machine = QtStateMachine.QStateMachine(self)

        self.state = {}
        self.state['idle'] = QtStateMachine.QState()
        self.state['idle_none'] = QtStateMachine.QState(self.state['idle'])
        self.state['idle_open'] = QtStateMachine.QState(self.state['idle'])
        self.state['idle_done'] = QtStateMachine.QState(self.state['idle'])
        self.state['idle_unchanged'] = QtStateMachine.QState(self.state['idle_done'])
        self.state['idle_updated'] = QtStateMachine.QState(self.state['idle_done'])
        self.state['idle_last'] = QtStateMachine.QHistoryState(self.state['idle'])
        self.state['idle'].setInitialState(self.state['idle_none'])
        self.state['idle_done'].setInitialState(self.state['idle_unchanged'])
        self.state['running'] = QtStateMachine.QState()

        self.machine.addState(self.state['idle'])
        self.machine.addState(self.state['running'])
        self.machine.setInitialState(self.state['idle'])

        state = self.state['idle']
        state.assignProperty(self.action['run'], 'visible', True)
        state.assignProperty(self.action['stop'], 'visible', False)
        state.assignProperty(self.action['open'], 'enabled', True)
        state.assignProperty(self.action['save'], 'enabled', True)

        state = self.state['idle_none']
        state.assignProperty(self.action['run'], 'enabled', False)
        state.assignProperty(self.action['save'], 'enabled', False)
        state.assignProperty(self.paramWidget, 'enabled', False)
        state.assignProperty(self.pane['param'], 'enabled', False)
        state.assignProperty(self.pane['list'], 'enabled', False)
        state.assignProperty(self.pane['list'].labelFoot, 'text', 'Nothing to show')
        state.assignProperty(self.pane['preview'], 'enabled', False)

        state = self.state['idle_open']
        state.assignProperty(self.action['run'], 'enabled', True)
        state.assignProperty(self.action['save'], 'enabled', False)
        state.assignProperty(self.paramWidget, 'enabled', True)
        state.assignProperty(self.pane['param'], 'enabled', True)
        state.assignProperty(self.pane['list'], 'enabled', False)
        state.assignProperty(self.pane['list'].labelFoot, 'text', 'Nothing to show')
        state.assignProperty(self.pane['preview'], 'enabled', False)

        state = self.state['idle_done']
        state.assignProperty(self.action['run'], 'enabled', True)
        state.assignProperty(self.action['save'], 'enabled', True)
        state.assignProperty(self.paramWidget, 'enabled', True)
        state.assignProperty(self.pane['param'], 'enabled', True)
        state.assignProperty(self.pane['list'], 'enabled', True)
        state.assignProperty(self.pane['list'].labelFoot, 'text', 'Double-click for preview')
        state.assignProperty(self.pane['preview'], 'enabled', True)
        def onEntry(event):
            self.splitter.setSizes([1, 1, -1])
        state.onEntry = onEntry

        state = self.state['running']
        state.assignProperty(self.action['run'], 'visible', False)
        state.assignProperty(self.action['stop'], 'visible', True)
        state.assignProperty(self.action['open'], 'enabled', False)
        state.assignProperty(self.action['save'], 'enabled', False)
        state.assignProperty(self.paramWidget, 'enabled', False)
        state.assignProperty(self.pane['param'], 'enabled', False)
        state.assignProperty(self.pane['list'], 'enabled', False)
        state.assignProperty(self.pane['preview'], 'enabled', False)

        state = self.state['idle_unchanged']
        def onEntry(event):
            self.pane['param'].flag = None
            self.pane['list'].flag = None
            self.pane['preview'].flag = None
            self.pane['param'].flagTip = None
            self.pane['list'].flagTip = None
            self.pane['preview'].flagTip = None
        state.onEntry = onEntry

        state = self.state['idle_updated']
        def onEntry(event):
            tip = ( 'Parameters have changed,\n' +
                    're-run analysis to update results.')
            self.pane['param'].flag = '*'
            self.pane['list'].flag = '*'
            self.pane['preview'].flag = '*'
            self.pane['param'].flagTip = tip
            self.pane['list'].flagTip = tip
            self.pane['preview'].flagTip = tip
        state.onEntry = onEntry

        transition = self.taggedTransition('OPEN')
        def onTransition(event):
            event = machine.TaggedEvent(event)
            file = event.kwargs['file']
            fileInfo = QtCore.QFileInfo(file)
            fileName = fileInfo.fileName()
            absolute = fileInfo.absoluteFilePath()
            self.setWindowTitle(self.title + ' - ' + fileName)
            self.line.file.setText(file)
            self.analysis.file = absolute
            self.folder.clear()
            self.pane['preview'].title = 'Preview'
            self.stack.setCurrentWidget(self.preview)
            self.preview.clear()
        transition.onTransition = onTransition
        transition.setTargetState(self.state['idle_open'])
        self.state['idle'].addTransition(transition)

        transition = self.taggedTransition('RUN')
        transition.setTargetState(self.state['running'])
        self.state['idle'].addTransition(transition)

        transition = self.taggedTransition('DONE')
        def onTransition(event):
            event = machine.TaggedEvent(event)
            self.folder.open(self.temp.name + '/')
            self.folder.setCurrentItem(self.folder.item(0))
            self.handlePreview(self.folder.item(0))
            msgBox = QtWidgets.QMessageBox(self)
            msgBox.setWindowTitle(self.windowTitle())
            msgBox.setIcon(QtWidgets.QMessageBox.Information)
            msgBox.setText('Analysis complete.')
            msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
            msgBox.setDefaultButton(QtWidgets.QMessageBox.Ok)
            self.msgShow(msgBox)
        transition.onTransition = onTransition
        transition.setTargetState(self.state['idle_unchanged'])
        self.state['running'].addTransition(transition)

        transition = self.taggedTransition('UPDATE')
        transition.setTargetState(self.state['idle_updated'])
        self.state['idle_done'].addTransition(transition)

        transition = self.taggedTransition('FAIL')
        def onTransition(event):
            event = machine.TaggedEvent(event)
            self.folder.open(self._temp.name + '/')
            self.folder.setCurrentItem(self.folder.item(0))
            self.handlePreview(self.folder.item(0))
            self.showException(*event.args)
        transition.onTransition = onTransition
        transition.setTargetState(self.state['idle_done'])
        self.state['running'].addTransition(transition)

        transition = self.taggedTransition('CANCEL')
        transition.setTargetState(self.state['idle_last'])
        self.state['running'].addTransition(transition)

        self.machine.start()

    def handleOpen(self, checked=False, fileName=None):
        """Called by toolbar action: open"""
        if fileName is None:
            fileName, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                self.title + ' - Open File',
                str(pathlib.Path.cwd()),
                'All Files (*) ;; Comma Separated Vector files (*.csv)')
        if len(fileName) == 0:
            return

        self.analysis = core.PartitionAnalysis(fileName)
        suffix = QtCore.QFileInfo(fileName).suffix()
        self.analysis.params.general.mega = (suffix == 'csv')

        self.setParamModel(self.analysis)
        self.postAction('OPEN', file=fileName)

    def handleSave(self):
        """Called by toolbar action: save"""
        path = pathlib.Path(self.analysis.file)

        # Dialog name filters with their file-selecting functions
        def fromSelection():
            return [pathlib.Path(item.file) for item in self.folder.selectedItems()]

        def fromFilter(filter):
            def _fromFilter():
                return list(pathlib.Path(self.temp.name).glob(filter))
            return _fromFilter

        nameFiltersWithSelectors = {
            'All files (*)': fromFilter('*.*'),
            'Selected files (*)': fromSelection,
            'Table files (*.tab)': fromFilter('*.tab'),
            'Spart files (*.spart)': fromFilter('*.spart'),
            'Partition files (*.txt)': fromFilter('*.txt'),
            'Tree files (*.tree)': fromFilter('*.tree'),
            'Vector Graphics (*.svg)': fromFilter('*.svg'),
            'Log files (*.log)': fromFilter('*.log'),
            }

        # Widget-based dialog, filters decide what files are saved
        dialog = QtWidgets.QFileDialog()
        dialog.setWindowTitle('Save with prefix')
        dialog.selectFile(str(path.stem))
        dialog.setDirectory(str(path.parent))
        dialog.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dialog.setAcceptMode(QtWidgets.QFileDialog.AcceptSave)
        dialog.setNameFilters(nameFiltersWithSelectors.keys())

        # Disable file selection
        class ProxyModel(QtCore.QIdentityProxyModel):
            def flags(self, index):
                flags = super().flags(index)
                if not self.sourceModel().isDir(index):
                    flags &= ~QtCore.Qt.ItemIsSelectable
                    flags &= ~QtCore.Qt.ItemIsEnabled
                return flags
        proxy = ProxyModel(dialog)
        dialog.setProxyModel(proxy)

        # All files will have this as prefix
        saveTo = ''
        if (dialog.exec()):
            saveTo = dialog.selectedFiles()[0]
            print('> Saving files to folder:',saveTo)
        if len(saveTo) == 0:
            return
        save = pathlib.Path(saveTo)

        # Select the files that will be copied
        filesFrom = nameFiltersWithSelectors[dialog.selectedNameFilter()]()

        filesMap = {
            file: save.with_name(save.name + '.' + file.name)
            for file in filesFrom
            }

        # Check existing files for possible overwrites
        existing = set(save.parent.glob('*'))
        overlap = existing & set(filesMap.values())

        if len(overlap) > 0:
            msgBox = QtWidgets.QMessageBox(self)
            msgBox.setWindowTitle(self.windowTitle())
            msgBox.setIcon(QtWidgets.QMessageBox.Question)
            msgBox.setText('Some files already exist. Overwrite?')
            msgBox.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
            msgBox.setDefaultButton(QtWidgets.QMessageBox.Yes)
            confirm = self.msgShow(msgBox)
            if confirm == QtWidgets.QMessageBox.No:
                return

        # Finally copy the files
        for file in filesMap.keys():
            print(file, '->', filesMap[file])
            shutil.copyfile(file, filesMap[file])

    def handleRun(self):
        """Called by toolbar action: run"""
        try:
            self._temp = tempfile.TemporaryDirectory(prefix='asap_')
            self.analysis.target = pathlib.Path(self._temp.name).as_posix()
        except Exception as exception:
            self.fail(exception)
            return

        def done(result):
            self.temp = self._temp
            self.analysis.results = result
            self.postAction('DONE', True)

        def fail(exception, trace):
            self.postAction('FAIL', exception, trace)

        self.launcher = threading.Process(work_run, self.analysis)
        self.launcher.done.connect(done)
        self.launcher.fail.connect(fail)
        # self.launcher.setLogger(logging.getLogger())
        self.launcher.start()
        self.postAction('RUN')

    def handleStop(self):
        """Called by cancel button"""
        msgBox = QtWidgets.QMessageBox(self)
        msgBox.setWindowTitle(self.windowTitle())
        msgBox.setIcon(QtWidgets.QMessageBox.Question)
        msgBox.setText('Cancel ongoing analysis?')
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No)
        msgBox.setDefaultButton(QtWidgets.QMessageBox.No)
        confirm = self.msgShow(msgBox)
        if confirm == QtWidgets.QMessageBox.Yes:
            self.launcher.quit()
            self.postAction('CANCEL')

    def handlePreview(self, item):
        """Called by file double-click"""
        if item is None:
            return
        try:
            path = pathlib.Path(item.file)
            self.pane['preview'].footer = path.name
            self.pane['preview'].title = 'Preview - ' + path.name
            if path.suffix == '.svg':
                self.stack.setCurrentWidget(self.graph)
                self.graphSvg.load(widgets.VectorPixmap.loadAndMap(str(path), self.colormap_graph[path.stem]))
                self.graphSvg.renderer().setAspectRatioMode(QtCore.Qt.KeepAspectRatio)
            else:
                self.stack.setCurrentWidget(self.preview)
                self.preview.clear()
                with open(item.file) as input:
                    for line in input:
                        self.preview.insertPlainText(line)
                format = QtGui.QTextBlockFormat()
                # format.setLineHeight(200, QtGui.QTextBlockFormat.ProportionalHeight)
                # format.setNonBreakableLines(True)
                format.setTopMargin(5)
                cursor = self.preview.textCursor()
                cursor.select(QtGui.QTextCursor.Document);
                cursor.mergeBlockFormat(format);
                self.preview.moveCursor(QtGui.QTextCursor.Start)
        except Exception as exception:
            self.fail(exception)
            return

    def showException(self, exception, trace):
        print(trace)
        msgBox = QtWidgets.QMessageBox(self)
        msgBox.setWindowTitle(self.title)
        msgBox.setIcon(QtWidgets.QMessageBox.Critical)
        msgBox.setText('An exception occured:')
        msgBox.setInformativeText(str(exception))
        msgBox.setDetailedText(trace)
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        self.msgShow(msgBox)
