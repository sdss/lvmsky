#!/usr/bin/env python
#
# ESO in kind team Innsbruck
# author:    Cesary Szyszka
# last edit: Stefan Kimeswenger
# 
#
#  This file is part of the SKYCORR software package.
#  Copyright (C) 2009-2013 European Southern Observatory
#
#  This programme is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This programme is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this programme. If not, see <http://www.gnu.org/licenses/>.
#
__version__ = '1.1.2'

# Used to guarantee to use at least Wx2.8
import wxversion
wxversion.ensureMinimal('2.8')

import matplotlib
import wx
import sys
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure
from matplotlib.ticker import NullFormatter
#from pylab import MultipleLocator

VERSION = 'SkyCorr GUI v%s' % __version__
global frame

class Model():
    def __init__(self):
        self.file_path = ''
        self.hdulist = None
        self.wavelength = []
        self.observed = []
        self.model = []
        self.residual = []

    def define_input(self, *args):
        import os.path
        if len(args) > 0:
            file_path = args[0]
        else:
            text = 'ERROR: no argument (FITS file) given in parameter list'\
                   '\n\n       EXIT'
            dlg = wx.MessageDialog( frame, text,
	                VERSION, wx.OK | wx.NO_DEFAULT)
            dlg.ShowModal() 
            sys.exit(-1)    
            
        if os.path.exists(file_path) and os.path.isfile(file_path) and \
           os.access(file_path, os.R_OK):
            self.file_path = file_path
        else:
            text = 'ERROR: FITS file \n"%s"\n'\
                   'given in parameter either does not exist or user "%s" ' \
                   'has no read permissions.\n\n'\
                   '       EXIT' % (file_path, os.getenv('LOGNAME'))
            dlg = wx.MessageDialog( frame, text,
	                VERSION, wx.OK | wx.NO_DEFAULT)
            dlg.ShowModal() 
            sys.exit(-1)

    def read_output_table(self):
        try:
            import pyfits
        except ImportError:
            import astropy.io.fits as pyfits

        # read in the data from the fits file and close it
        self.hdulist = pyfits.open(self.file_path)
        #self.hdulist.info()
        _data = self.hdulist[1].data
        _colum_names = self.hdulist[1].columns.names
        self.hdulist.close()

        # parse the data into correct fields
        for line in _data:
            self.wavelength.append(line[_colum_names.index('lambda')])
            self.observed.append(line[_colum_names.index('flux')])
            self.model.append(line[_colum_names.index('mflux')])
            self.residual.append(line[_colum_names.index('scflux')])


#class SkyCorrFrame(wx.Frame):
#    def __init__(self, *args, **kwargs):

class SkyCorrFrame(wx.Frame):
    def __init__(self, *args, **kwargs):
        wx.Frame.__init__(self, *args, **kwargs)
        self.SetBackgroundColour(wx.NamedColour("WHITE"))
        self.figure = Figure()
        self.canvas = FigureCanvas(self, -1, self.figure)
        self.lines = []

        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        tw, th = self.toolbar.GetSizeTuple()
        fw, fh = self.canvas.GetSizeTuple()
        self.toolbar.SetSize(wx.Size(fw, th))

        # removed here as we have no data yet
        #self.draw()

        # exit button
        self.button = wx.Button(self, -1, "Exit")
        self.Bind(wx.EVT_BUTTON, self.OnExit, self.button)

        # text for measuring with the mouse
        self.position_text = wx.StaticText(
            self,-1, "                                      "
                     "                                    \n", 
            (40,110),(-1,-1), wx.ALIGN_LEFT)

        # Bottom sizer (Horizontal)
        self.bottom_sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.bottom_sizer.Add(self.toolbar, 3, wx.EXPAND)
        self.bottom_sizer.AddSpacer(30)
        self.bottom_sizer.Add(self.position_text, 0, border=3, 
            flag=wx.ALIGN_LEFT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        self.bottom_sizer.AddSpacer(30)
        self.bottom_sizer.Add(self.button, 1, wx.EXPAND)
        self.bottom_sizer.AddSpacer(10)

        # Vertical sizer
        self.vertical_sizer = wx.BoxSizer(wx.VERTICAL)
        self.vertical_sizer.Add(self.canvas, 1, wx.EXPAND)
        self.vertical_sizer.Add(self.bottom_sizer, 0, wx.EXPAND)

        # set sizer
        self.SetSizer(self.vertical_sizer)

        # setup menu bar
        file_menu = wx.Menu()
        menuAbout = file_menu.Append(wx.ID_ABOUT, '&About', \
                    'Information about this program.')
        file_menu.AppendSeparator()
        menuExit = file_menu.Append(wx.ID_EXIT, 'E&xit', \
            'Terminate the program')
        menu_bar = wx.MenuBar()
        menu_bar.Append(file_menu, '&File')

        # Actually add the menu bar
        self.SetMenuBar(menu_bar)

        # bind handlers with menu items
        self.Bind(wx.EVT_MENU, self.OnAbout, menuAbout)
        self.Bind(wx.EVT_MENU, self.OnExit, menuExit)

    def OnAbout(self, event):
        # Show a short message about this program. 
        dialog = wx.MessageDialog(self,  
         'This GUI is the plot tool for the SM-02 SkyCorr program.\n\n' +
         'M. Barden, W. Kausch, S. Kimeswenger (head),\n' +
         'A.M. Jones, S. Noll, C. Szyszka (GUI plot tool author)\n' +
         'Austrian ESO In-Kind Team Innsbruck\n\n' +
         'Institute for Astro and Particle Physics\n' +
         'University of Innsbruck, Austria\n\n' +
         '(c) ESO 2013', VERSION, wx.OK)
        dialog.ShowModal()
        dialog.Destroy()

    def OnExit(self, event):
        self.Close(True)
        sys.exit(0)   

    def sizeHandler(self, *args, **kwargs):
        self.canvas.SetSize(self.GetSize())
        self.toolbar.update()

    def draw(self):
        if not hasattr(self, 'subplot1'):

            # make shared x-axis
            self.figure.subplots_adjust(hspace=0.10)
            self.subplot1 = self.figure.add_subplot(211)
            self.subplot2 = self.figure.add_subplot(212)
            #self.subplot1.xaxis.set_major_formatter(NullFormatter())

        # Labels
        self.subplot1.set_title("SkyCorr results")
        self.subplot1.set_ylabel("Flux")
        self.subplot2.set_ylabel("Flux (observation - calculation)")
        self.subplot2.set_xlabel('Wavelength [microns]')

        self.lines += self.subplot1.plot(model.wavelength, model.observed, \
            color='black', linewidth=1,label='Input Science')
        self.lines += self.subplot1.plot(model.wavelength, model.model, \
            color='red', linewidth=1, label='Calculated Model')
        self.lines += self.subplot2.plot(model.wavelength, model.residual, 
            color='blue', linewidth=1, label='Sky-subtracted science')
        self.subplot1.legend(loc=2)
        self.subplot2.legend(loc=2)
        # add events for measuring in the plot
        self.canvas.mpl_connect('axes_enter_event'   , self.enter_graph)
        self.canvas.mpl_connect('axes_leave_event'   , self.exit_graph)


    def enter_graph(self,event):
        # when entering the graph the motion is monitored to display the data
        self.cid_move = self.canvas.mpl_connect('motion_notify_event', 
                                                 self.move_graph)

    def exit_graph(self,event):
        # on exit the data value refresh is turned off
        self.canvas.mpl_disconnect(self.cid_move)
        self.position_text.SetLabel("")    

    def move_graph(self,event):
        # refresh the data display itself
        if event.xdata != None and event.ydata != None:
            #print "measure:", event.xdata, event.ydata
            #if (GUI_Global.Y[GUI_Data.pixel_of_data(event.xdata)] != None and
            #    GUI_Data.pixel_of_data(event.xdata) != None):
            text = (('wavelength: %.5g  y: %.5g\n') % \
                    (event.xdata,event.ydata))
            #    (event.xdata,event.ydata,GUI_Data.pixel_of_data(event.xdata)+1,
            #    GUI_Global.Y[GUI_Data.pixel_of_data(event.xdata)]))
            #else:
            #    text = 'wavelength: %g     y: %g' % (event.xdata,event.ydata)
            self.position_text.SetLabel(text)    


class App(wx.App):

    def OnInit(self):
        global frame
        # Create the main window and insert the custom frame
        frame = SkyCorrFrame(parent=None, title=VERSION, \
            size=(1024, 724), pos=(40,40), style=wx.MINIMIZE_BOX | wx.CAPTION)
        return True


# start the application
if __name__ == "__main__":

    # initialize model object
    model = Model()

    # first start the mini app - otherwhise error messages with popups/dialogs
    # are not possible
    app = App(False)

    # Define and read input file after the main frame exists
    if len(sys.argv) > 1:
        model.define_input(sys.argv[1])
    else:
        # No input parameter given!
        model.define_input()

    # read table
    model.read_output_table()

    # draw and show after the data was read
    frame.Show(True)
    SkyCorrFrame.draw(frame)

    app.MainLoop()
