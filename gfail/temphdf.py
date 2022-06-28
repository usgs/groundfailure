#!/usr/bin/env python3

import os
import tables
import numpy as np

from mapio.shake import ShakeGrid


class TempHdf(object):
    def __init__(self, grid2dfile, filename, name=None):
        """
        Convert grid2d file into a temporary hdf5 file for reducing memory
        load.

        Args:
            grid2dfile: grid2d file object to save
            filename (str): Path to where file should be saved (recommended
                it be a temporary dir).
            name (str): Name of layer, if None, will use filename minus the
                extension, or if a multihazard grid2d object, each layer will
                have its own name.
        """
        filename1, file_ext = os.path.splitext(filename)
        if file_ext != ".hdf5":
            filename = filename1 + ".hdf5"
            print("Changed extension from %s to .hdf5" % (file_ext,))
        filters = tables.Filters(complevel=5, complib="blosc")
        with tables.open_file(filename, mode="w") as self.tempfile:
            self.gdict = grid2dfile.getGeoDict()
            if type(grid2dfile) == ShakeGrid:
                for layer in grid2dfile.getLayerNames():
                    filldat = grid2dfile.getLayer(layer).getData()
                    self.tempfile.create_carray(
                        self.tempfile.root, name=layer, obj=filldat, filters=filters
                    )
                self.shakedict = grid2dfile.getShakeDict()
                self.edict = grid2dfile.getEventDict()
            else:
                if name is None:
                    name = os.path.basename(filename1)
                filldat = grid2dfile.getData()
                self.tempfile.create_carray(
                    self.tempfile.root, name=name, obj=filldat, filters=filters
                )
            self.filename = os.path.abspath(filename)

    def getFilepath(self):
        """
        Return file path.
        """
        return self.filename

    def getGeoDict(self):
        """
        Return geodictionary.
        """
        return self.gdict

    def getShakeDict(self):
        """
        Return shake dictionary if it exists.
        """
        try:
            return self.shakedict
        except Exception as e:
            print(e)
            print("no shake dictionary found")
            return None

    def getEventDict(self):
        """
        Return event dictionary if it exists.
        """
        try:
            return self.edict
        except Exception as e:
            print(e)
            print("no event dictionary found")
            return None

    def getSlice(
        self, rowstart=None, rowend=None, colstart=None, colend=None, name=None
    ):
        """
        Return specified slice of data.

        Args:
            rowstart (int, None): Starting row index (inclusive), if None, will
                start at 0.
            rowend (int, None): Ending row index (exclusive), if None, will
                end at last row.
            colstart (int, None): Starting column index (inclusive), if None,
                will start at 0.
            colend (int, None): Ending column index (exclusive), if None, will
                end at last row.
            name (str): Name of layer/child name to return.

        Returns:
            numpy array of data.
        """
        if name is None:
            name, ext = os.path.splitext(os.path.basename(self.getFilepath()))
        if rowstart is None:
            rowstart = ""
        else:
            rowstart = int(rowstart)
        if rowend is None:
            rowend = ""
        else:
            rowend = int(rowend)
        if colstart is None:
            colstart = ""
        else:
            colstart = int(colstart)
        if colend is None:
            colend = ""
        else:
            colend = int(colend)

        indstr = "%s:%s, %s:%s" % (rowstart, rowend, colstart, colend)
        # So end index will actually be captured:
        indstr = indstr.replace("-1", "")
        with tables.open_file(self.filename, mode="r") as file1:
            try:
                dataslice = eval("file1.root.%s[%s]" % (name, indstr))
                return dataslice
            except Exception as e:
                raise Exception(e)
        return

    def getSliceDiv(self, rowmax=None, colmax=None):
        """
        Determine how to slice the arrays.

        Args:
            rowmax (int): Maximum number of rows in each slice; default None
                uses entire row.
            colmax (int): Maximum number of columns in each slice; default
                None uses entire column.

        Returns:
            tuple: rowstarts, rowends, colstarts, colends.
        """
        numrows = self.gdict.ny
        numcols = self.gdict.nx
        if rowmax is None or rowmax > numrows:
            rowmax = numrows
        if colmax is None or colmax > numcols:
            colmax = numcols
        numrowslice, rmrow = divmod(numrows, rowmax)
        numcolslice, rmcol = divmod(numcols, colmax)
        rowst = np.arange(0, numrowslice * rowmax, rowmax)
        rowen = np.arange(rowmax, (numrowslice + 1) * rowmax, rowmax)
        if rmrow > 0:
            rowst = np.hstack([rowst, numrowslice * rowmax])
            rowen = np.hstack([rowen, None])
        else:
            rowen = np.hstack([rowen[:-1], None])
        colst = np.arange(0, numcolslice * colmax, colmax)
        colen = np.arange(colmax, (numcolslice + 1) * colmax, colmax)
        if rmcol > 0:
            colst = np.hstack([colst, numcolslice * colmax])
            colen = np.hstack([colen, None])
        else:
            colen = np.hstack([colen[:-1], None])
        rowstarts = np.tile(rowst, len(colst))
        colstarts = np.repeat(colst, len(rowst))
        rowends = np.tile(rowen, len(colen))
        colends = np.repeat(colen, len(rowen))
        return rowstarts, rowends, colstarts, colends
