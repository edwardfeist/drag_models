import numpy as np
import plotly.graph_objects as go
import argparse

def plot_model(fig, mname, model, Mach, Re, T_p, T_inf, gam, Kn):
    
    name=f"{mname}: Mach={Mach}"
    Cd = model(Mach, Re, T_p, T_inf, gam, Kn)
    
    fig.add_trace(go.Scatter(x=Re,
                           y=Cd,
                           mode="lines",
                           name=name,
                           meta=[name],
                           hovertemplate="<b>%{meta[0]}</b>" + 
                                           "<br>Re: %{x:.3e}" +
                                           "<br>Cd: %{y:.3e}" +
                                           "<extra></extra>"))
    
def plot_subfunction(fig, mname, model, Mach, Re, T_p, T_inf, gam, Kn):
    
    name=f"{mname}: Re={Re}"
    if not isinstance(Mach, np.ndarray):
        Mach = np.asarray([Mach])
        
    Cd=np.zeros(np.shape(Mach))
    for rei in np.arange(len(Mach)):
        Cd[rei] = model(Mach[rei])
    
    fig.add_trace(go.Scatter(x=Mach,
                           y=Cd,
                           mode="lines",
                           name=name,
                           meta=[name],
                           hovertemplate="<b>%{meta[0]}</b>" + 
                                           "<br>Re: %{x:.3e}" +
                                           "<br>Cd: %{y:.3e}" +
                                           "<extra></extra>"))
    
def setlayout(fig):
    fig.update_layout(xaxis_type="log",
                       yaxis_type="log",
                       font_size=18,
                       font=dict(family="Times New Roman", color="black"),
                       plot_bgcolor="white",
                       showlegend=True,
                       )
    fig.update_xaxes(title=r"$Re_p$",
                      showline=True,
                      mirror=True,
                      linewidth=2,
                      linecolor="black",
                      showgrid=True,
                      gridcolor="lightgrey",
                      )
    fig.update_yaxes(title="$C_d$",
                      showline=True,
                      mirror=True,
                      linewidth=2,
                      linecolor="black",
                      showgrid=True,
                      gridcolor="lightgrey",
                      zeroline=True,
                      zerolinecolor="lightgrey",
                      ) 
    
def list_of_paths(arg):
    return arg.split(',')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot drag models. Arguments are comma-separated lists."
        )
    
    parser.add_argument("--plot", type=list_of_paths, default=None)
    parser.add_argument("--mach", type=list_of_paths, default=None)
    
    args=parser.parse_args()
    
    if args.plot is not None:
        
        plotlist = [x.lower() for x in args.plot]
        
        print(f"Requested list of plots: {plotlist}")
            
        Re = np.logspace(0, 6, 1000)
        T_p = 1000
        T_inf = 1000
        gam = 1.4
        Kn = 0.00001
        machlist = [0.24, 0.39, 0.79, 1.18, 1.57, 2.36, 3.15, 4.72]
        
        if args.mach is not None:
            machlist = [float(x) for x in args.mach]
        
        if "henderson" in plotlist:
            from sphere.henderson import Henderson_Cd
            figH=go.Figure()
            
            for m in machlist:
                plot_model(figH, "Henderson", Henderson_Cd, m, Re, T_p, T_inf, gam, Kn)
            
            setlayout(figH)  
            figH.show()
            
        if "loth2008" in plotlist:
            from sphere.loth2008 import Loth_2008_Cd
            figL=go.Figure()
            
            for m in machlist:
                plot_model(figL, "Loth 2008", Loth_2008_Cd, m, Re, T_p, T_inf, gam, Kn)
            
            setlayout(figL)  
            figL.show()
            
        if "loth2021" in plotlist:
            from sphere.loth2021 import Loth_2021_Cd
            figLL=go.Figure()
            
            for m in machlist:
                plot_model(figLL, "Loth 2021", Loth_2021_Cd, m, Re, T_p, T_inf, gam, Kn)
            
            setlayout(figLL)  
            figLL.show()
            
        if "test_2021_subfunctions" in plotlist:
            from sphere.loth2021 import C21, G21, H21, J21
            figLT=go.Figure()
            
            Mach = np.linspace(0.1,10,1000)
            plot_subfunction(figLT, "Cm", C21, Mach, 45, T_p, T_inf, gam, Kn)
            plot_subfunction(figLT, "Gm", G21, Mach, 45, T_p, T_inf, gam, Kn)
            plot_subfunction(figLT, "Hm", H21, Mach, 45, T_p, T_inf, gam, Kn)
            plot_subfunction(figLT, "Jm", J21, Mach, 45, T_p, T_inf, gam, Kn)
            
            setlayout(figLT) 
            figLT.update_xaxes(title=r"$Mach$",
                              )
            figLT.update_yaxes(title="$Function Value$",
                              ) 
            figLT.show()