# 
# Maillages autour de formes simples
#
import Generator as G
__version__ = G.__version__

def square((x0,y0), (x1,y1), H, (ni,nj,nk)):
    """Return the 8 arrays for a cube"""
    out = []
    
    a = G.cart((x0-H,y0,0.), (H*1./(nk-1), 
                              (y1-y0)*1./(nj-1), 1.), (nk,nj,1))
    out.append(a)

    a = G.cart( (x0,y1,0.), ((x1-x0)*1./(ni-1), 
                             H*1./(nk-1), 1.), (ni,nk,1))
    out.append(a)

    a = G.cart( (x1,y0,0.), (H*1./(nk-1), 
                             (y1-y0)*1./(nj-1), 1.), (nk,nj,1))
    out.append(a)
    
    a = G.cart( (x0,y0-H,0.), ((x1-x0)*1./(ni-1), 
                               H*1./(nk-1), 1.), (ni,nk,1))
    out.append(a)

    a = G.cart( (x0-H,y0-H,0.), (H*1./(nk-1), 
                                 H*1./(nk-1), 1.), (nk,nk,1))
    out.append(a)

    a = G.cart( (x0-H,y1,0.), (H*1./(nk-1), 
                               H*1./(nk-1), 1.), (nk,nk,1))
    out.append(a)

    a = G.cart( (x1,y1,0.), (H*1./(nk-1), 
                             H*1./(nk-1), 1.), (nk,nk,1))
    out.append(a)

    a = G.cart( (x1,y0-H,0.), (H*1./(nk-1), 
                               H*1./(nk-1), 1.), (nk,nk,1))
    out.append(a)

    return out
