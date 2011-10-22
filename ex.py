from openeye.oechem import *
import stats

from matplotlib.figure import Figure
from matplotlib.patches import Polygon
from matplotlib.backends.backend_agg import FigureCanvasAgg
import matplotlib.numerix as nx

def Ellipse((x,y), (rx, ry), resolution=20, orientation=0, **kwargs):
    theta = 2*nx.pi/resolution*nx.arange(resolution) + orientation
    xs = x + rx * nx.cos(theta)
    ys = y + ry * nx.sin(theta)
    return Polygon(zip(xs, ys), **kwargs)

# Read up to 'limit' records that have XLogP values
# Return the lists of:
#  identifiers, molecular weights, XLogP values
def read_data(ifs, limit = None):
    cids = []
    weights = []
    xlogps = []
    for i, mol in enumerate(ifs.GetOEGraphMols()):
        # Some of the compounds don't have an XLOGP value
        # Skip those molecules
        if not OEHasSDData(mol, "PUBCHEM_CACTVS_XLOGP"):
            continue
        cid = OEGetSDData(mol, "PUBCHEM_COMPOUND_CID")
        weight = OEGetSDData(mol, "PUBCHEM_OPENEYE_MW")
        xlogp = OEGetSDData(mol, "PUBCHEM_CACTVS_XLOGP")
        if (cid == "" or weight == "" or xlogp == ""):
            raise AssertionError( (cid, weight, xlogp) )

        cids.append(cid)
        weights.append(float(weight))
        xlogps.append(float(xlogp))

        if limit is not None and len(cids) >= limit:
            break
        
    return cids, weights, xlogps

def calculate_ellipse_data(xdata, ydata):
    xcenter = stats.lmean(xdata)
    xradius = stats.lstdev(xdata)
    ycenter = stats.lmean(ydata)
    yradius = stats.lstdev(ydata)
    return (xcenter, ycenter), (xradius, yradius)

def main():
    filename = "/Users/dalke/databases/compounds_500001_510000.sdf.gz"
    ifs = oemolistream(filename)

    # The figure will be 3 inches by 3 inches
    # Ths size is important because the text is defined relative to
    # inches and not pixels.  In a smaller image the text is more
    # cramped and likely to overlap.  In a larger image the text is
    # not big enough.  This works well for my plot.
    fig = Figure(figsize=(4,4))
    
    ax = fig.add_subplot(111)

    cids, weights, xlogps = read_data(ifs, 100)
    ax.scatter(weights, xlogps)
    center, radii = calculate_ellipse_data(weights, xlogps)
    ax.add_patch(Ellipse(center, radii, fill=0, edgecolor="blue"))
    
    cids, weights, xlogps = read_data(ifs, 100)
    ax.scatter(weights, xlogps, marker = "^", color="red")
    center, radii = calculate_ellipse_data(weights, xlogps)
    ax.add_patch(Ellipse(center, radii, fill=0, edgecolor="red"))
    
    ax.set_xlabel("Atomic weight")
    ax.set_ylabel("CACTVS XLogP")

        # Make the PNG
    canvas = FigureCanvasAgg(fig)
    # The size * the dpi gives the final image size
    #   a4"x4" image * 80 dpi ==> 320x320 pixel image
    canvas.print_figure("mw_v_xlogp_ellipses.png", dpi=80)

if __name__ == "__main__":
    OEThrow.SetLevel(OEErrorLevel_Error)
    main()