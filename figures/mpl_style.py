"""
Color constants and matplotlib style definitions
"""

import matplotlib as mpl
import cmcrameri.cm as cmc
import cartopy.crs as ccrs


CMAP = cmc.batlow
CMAP_an = cmc.vik
CMAP_discr = cmc.batlowW


mpl.rcParams["legend.frameon"] = False
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = [
    "Tahoma",
    "DejaVu Sans",
    "Lucida Grande",
    "Verdana",
]

mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.grid"] = True

mpl.rcParams["savefig.dpi"] = 300
mpl.rcParams["savefig.bbox"] = "tight"
mpl.rcParams["savefig.transparent"] = True

def plot_cities_expats(ax, color, symbol_size):
    
    """add cities on a plot"""
    
    # set cities coordinates
    trento = [46.0667, 11.1167] #lat, lon
    bolzano = [46.4981, 11.3548]
    Munich = [48.137154, 11.576124]
    Innsbruck = [47.259659, 11.400375]
    Salzburg = [47.811195, 13.033229]
    Stuttgard = [48.783333, 9.183333]
    Zurich = [47.373878, 8.5450984]
    Ulm = [48.4, 9.98]
    Ravensburg = [47.7, 9.61]
    Regensburg = [49.01, 12.11]
    Nurnberg = [49.5, 11.08]
    Jena = [50.92, 11.59]
    Frankfurt = [50.174, 8.62]
    Freiburg = [48.059, 8.05]
    Padova = [45.4064, 11.8768]
    Venice = [45.4408, 12.3155]
    Milan = [45.4642, 9.1900]
    Bolzano = [46.4981, 11.3548]
    Bologna = [44.4949, 11.3426]
    Trento = [46.0667, 11.1167]
    Ferrara = [44.8353, 11.6198]
    Ravenna = [44.4184, 12.2035]

    # ass swiss cities
    Bern = [46.9481, 7.4474]
    Lausanne = [46.5197, 6.6323]
    Zurich = [47.373878, 8.5450984]

    
    # Plot the points
    ax.scatter(Ulm[1], Ulm[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Munich[1], Munich[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Stuttgard[1], Stuttgard[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Ravensburg[1], Ravensburg[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Regensburg[1], Regensburg[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Nurnberg[1], Nurnberg[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Jena[1], Jena[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Frankfurt[1], Frankfurt[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Freiburg[1], Freiburg[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Padova[1], Padova[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Venice[1], Venice[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Milan[1], Milan[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Bolzano[1], Bolzano[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Bologna[1], Bologna[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Trento[1], Trento[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Ferrara[1], Ferrara[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Ravenna[1], Ravenna[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())    
    ax.scatter(Bern[1], Bern[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Lausanne[1], Lausanne[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Zurich[1], Zurich[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    
    # Plot the names next to the points, adjusted for lower right positioning
    ax.text(Ulm[1] + 0.02, Ulm[0] - 0.02, 'Ulm', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Munich[1] + 0.02, Munich[0] - 0.02, 'Munich', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Stuttgard[1] + 0.02, Stuttgard[0] - 0.02, 'Stuttgard', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Ravensburg[1] + 0.02, Ravensburg[0] - 0.02, 'Ravensburg', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Regensburg[1] + 0.02, Regensburg[0] - 0.02, 'Regensburg', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Nurnberg[1] + 0.02, Nurnberg[0] - 0.02, 'Nurnberg', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Jena[1] + 0.02, Jena[0] - 0.02, 'Jena', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Frankfurt[1] + 0.02, Frankfurt[0] - 0.02, 'Frankfurt', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Freiburg[1] + 0.02, Freiburg[0] - 0.02, 'Freiburg', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Padova[1] + 0.02, Padova[0] - 0.02, 'Padova', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Venice[1] + 0.02, Venice[0] - 0.02, 'Venice', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Milan[1] + 0.02, Milan[0] - 0.02, 'Milan', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Bolzano[1] + 0.02, Bolzano[0] - 0.02, 'Bolzano', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Bologna[1] + 0.02, Bologna[0] - 0.02, 'Bologna', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Trento[1] + 0.02, Trento[0] - 0.02, 'Trento', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Ferrara[1] + 0.02, Ferrara[0] - 0.02, 'Ferrara', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Ravenna[1] + 0.02, Ravenna[0] - 0.02, 'Ravenna', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')    
    ax.text(Bern[1] + 0.02, Bern[0] - 0.02, 'Bern', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Lausanne[1] + 0.02, Lausanne[0] - 0.02, 'Lausanne', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Zurich[1] + 0.02, Zurich[0] - 0.02, 'Zurich', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')   
    
    return

    
def plot_cities_IT(ax, color, symbol_size):
    
    """add cities on a plot"""
    
    # set cities coordinates
    trento = [46.0667, 11.1167] #lat, lon
    bolzano = [46.4981, 11.3548]
    Padova = [45.4064, 11.8768]
    Venice = [45.4408, 12.3155]
    Milan = [45.4642, 9.1900]
    Bolzano = [46.4981, 11.3548]
    Bologna = [44.4949, 11.3426]
    Trento = [46.0667, 11.1167]
    Turin = [45.0703, 7.6869]
    Genoa = [44.4056, 8.9463]
    Florence = [43.7696, 11.2558]

    # Plot the points
    ax.scatter(Padova[1], Padova[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    #ax.scatter(Venice[1], Venice[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Milan[1], Milan[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Bolzano[1], Bolzano[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Bologna[1], Bologna[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Trento[1], Trento[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Turin[1], Turin[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Genoa[1], Genoa[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(Florence[1], Florence[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    
    # Plot the names next to the points, adjusted for lower right positioning
    ax.text(Padova[1] + 0.02, Padova[0] - 0.02, 'Padova', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    #ax.text(Venice[1] + 0.02, Venice[0] - 0.02, 'Venice', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Milan[1] + 0.02, Milan[0] - 0.02, 'Milan', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Bolzano[1] + 0.02, Bolzano[0] - 0.02, 'Bolzano', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Bologna[1] + 0.02, Bologna[0] - 0.02, 'Bologna', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Trento[1] + 0.02, Trento[0] - 0.02, 'Trento', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Turin[1] + 0.02, Turin[0] - 0.02, 'Turin', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Genoa[1] + 0.02, Genoa[0] - 0.02, 'Genoa', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Florence[1] + 0.02, Florence[0] - 0.02, 'Florence', color=color, transform=ccrs.PlateCarree(), ha='left', va='top') 
        
    return

    
    
def plot_local_dfg(ax, color, symbol_size):
    
    
    #add instrument positions 
    Penegal = [46.43921, 11.2155]
    Tarmeno = [46.34054, 11.2545]
    Vilpiano = [46.55285, 11.20195]
    Sarntal = [46.6412136, 11.3541082]
    Cles_Malgolo = [46.38098, 11.08136]
    trento = [46.0667, 11.1167] #lat, lon
    bolzano = [46.4981, 11.3548]  
    Soprabolzano = [46.5222, 11.3871]
    Corno_renon = [46.615, 11.460833]
    
     
    ax.scatter(trento[1], trento[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    ax.scatter(bolzano[1], bolzano[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    #ax.scatter(Penegal[1], Penegal[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())
    #ax.scatter(Cles_Malgolo[1], Cles_Malgolo[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())                        
    #ax.scatter(Tarmeno[1], Tarmeno[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())            
    ax.scatter(Vilpiano[1], Vilpiano[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())            
    ax.scatter(Sarntal[1], Sarntal[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())            
    ax.scatter(Soprabolzano[1], Soprabolzano[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())            
    ax.scatter(Corno_renon[1], Corno_renon[0], marker='x', color=color, s=symbol_size, transform=ccrs.PlateCarree())

    ax.text(trento[1] + 0.02, trento[0] - 0.02, 'Trento', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(bolzano[1] + 0.02, bolzano[0] - 0.02, 'Bolzano', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    #ax.text(Penegal[1] + 0.02, Penegal[0] - 0.02, 'Penegal', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    #ax.text(Cles_Malgolo[1] + 0.02, Cles_Malgolo[0] - 0.02, 'Cles_Malgolo', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    #ax.text(Tarmeno[1] + 0.02, Tarmeno[0] - 0.02, 'Tarmeno', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Vilpiano[1] + 0.02, Vilpiano[0] - 0.02, 'Vilpiano', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Sarntal[1] + 0.02, Sarntal[0] + 0.02, 'Sarntal', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Soprabolzano[1] + 0.02, Soprabolzano[0] - 0.02, 'Soprabolzano', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')
    ax.text(Corno_renon[1] + 0.02, Corno_renon[0] - 0.02, 'Corno del Renon', color=color, transform=ccrs.PlateCarree(), ha='left', va='top')

    return

