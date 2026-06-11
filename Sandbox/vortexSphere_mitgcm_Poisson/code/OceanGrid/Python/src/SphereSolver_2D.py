#%%
#COMMENT: Libraries imported
import csv, os ,sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix

from scipy.sparse.linalg import spsolve , bicgstab ,LinearOperator, spilu , gmres,cg

from scipy.sparse import csr_matrix

from src.config import (nLat, nPhi,Lat, sinLat0, dLat, GRAVITY ,D, sinLat, OMEGA,R, ScaleEquation, neighborhandling, Ld0, includeFreeSurface, OUTPUT_PATH)
import plotly.express as px
import plotly.graph_objects as go

from scipy.sparse import coo_matrix

from tqdm import tqdm

#%%
#ADDED CLASS: SPHERE SOLVER 2D 
class SphereSolver2D:
    def __init__(self, Ocean_Grid=None):
        self.d = Ocean_Grid.d
        self.Ocean_Grid = Ocean_Grid
        self.nLat = nLat
        self.nPhi = nPhi
        self.R = R
        self.OMEGA = OMEGA
        self.GRAVITY = GRAVITY
        self.D = D
        self.Ld0 = Ld0
        self.dLat = dLat
        self.Lat = Lat

        self.sinLat =sinLat
        self.sinLat0 =sinLat0
        self.CORIOLIS = 2.0*OMEGA*sinLat

    def Psi_build_matrix(self):
        print('BUILDING COEFFICIENT MATRIX FOR ψ')
        n_elements = self.nLat * self.nPhi
        A = lil_matrix((n_elements, n_elements))

        unscaled_printed = False 
        scaled_printed= False
        ocean_count=0
        count=0
        
        for iLat in range(self.nLat): 
            Lat = self.Lat[iLat]
            cosLat = np.cos(Lat)
            sinLat = np.sin(Lat)
            LatPHalf = Lat + self.dLat/2
            LatMHalf = Lat - self.dLat/2
            cosLatPHalf = np.cos(LatPHalf)
            cosLatMHalf = np.cos(LatMHalf)

            #ROCKET!!: Neumann BC           
            if iLat==0: 
                cosLatMHalf=0.0
            if iLat == self.nLat-1: 
                cosLatPHalf=0.0

            for iLon in range(self.nPhi):
                idx = iLat * self.nPhi + iLon

                #!Same as Fortran
                # if self.d[iLat,iLon] > 0:
                #     d = self.d[iLat,iLon]
                # else:
                #     d=0.0
               
                # #COMMENT: CENTER 
                # coef_center = d * (- 2.0  - (cosLatMHalf + cosLatPHalf) * cosLat ) - ((2.0 * OMEGA * sinLat)**2 /GRAVITY)* (R**2*self.dLat**2 *cosLat**2)

                # #COMMENT: SOUTH
                # if iLat > 0:
                #     coef_south  = cosLatMHalf  * cosLat *  d

                # #COMMENT: NORTH
                # if iLat < self.nLat-1:
                #     coef_north  = cosLatPHalf  * cosLat *  d

                # #COMMENT: WEST
                # coef_west   = 1.0 *  d
                    
                # #COMMENT: EAST
                # coef_east   = 1.0 *  d



                #!Neighbor handling
                #  # COMMENT: SKIP LAND CELL and continue the nested loop
                if self.d[iLat,iLon] <= 0:
                    coef_center = 1.0
                    A[idx, idx] = coef_center
                    continue

                #enters this chunk if the condition above is not satisfied meaning depth>0:
                #COMMENT: CENTER 
                coef_center = self.d[iLat,iLon]* (- 2.0  - (cosLatMHalf + cosLatPHalf) * cosLat ) - ((2.0 * OMEGA * sinLat)**2 /GRAVITY)* (R**2*self.dLat**2 *cosLat**2)

                coef_south  = 0.0
                coef_north  = 0.0
                coef_west   = 0.0
                coef_east   = 0.0

                #COMMENT: SOUTH; 2 CONDITIONS ON DEPTH AND ON SOUTH POLE 
                if iLat > 0 and self.d[iLat-1,iLon] > 0:
                    coef_south  = cosLatMHalf  * cosLat *  self.d[iLat,iLon]

                #COMMENT: NORTH ; 2 CONDITIONS ON DEPTH AND ON NORTH POLE
                if iLat < self.nLat-1 and self.d[iLat+1,iLon] > 0:
                    coef_north  = cosLatPHalf  * cosLat *  self.d[iLat,iLon]

                #COMMENT: WEST
                left_iLon = iLon-1 if iLon > 0 else self.nPhi-1
                if self.d[iLat,left_iLon]>0:
                    coef_west   = 1.0 *  self.d[iLat,iLon]
                    
                #COMMENT: EAST
                right_iLon = iLon+1 if iLon < self.nPhi-1 else 0
                if self.d[iLat,right_iLon]>0:
                    coef_east   = 1.0 *  self.d[iLat,iLon]
                    
                # #? Filling the matrix
                A[idx, idx] = coef_center

                # # #!ROCKET!! condition from fortran
                if iLat > 0:
                    A[idx , idx  - self.nPhi] = coef_south
                if iLat < self.nLat - 1:
                    A[idx , idx + self.nPhi  ] = coef_north

                # # #?wrapping around the longitude for the west-east
                left_idx = idx - 1; right_idx = idx + 1
                if iLon == 0:
                    left_idx = idx + (self.nPhi - 1)  
                elif iLon == self.nPhi - 1:
                    right_idx = idx - (self.nPhi - 1) 
                A[idx, left_idx]  = coef_west
                A[idx, right_idx] = coef_east 

                # - (R**2*self.dLat**2 *cosLat) * ((2*OMEGA*sinLat)**2/(GRAVITY*self.D))
                #+ (1.0/self.Ld0)**2 -2/(R**2*self.dLat**2 *cosLat**2) 
                # -2.0  - (cosLatMHalf +cosLatPHalf)* cosLat - ( (2.0*OMEGA*self.sinLat0)**2 * (self.R**2 * cosLat**2 * self.dLat**2) ) / ( GRAVITY * self.d[iLat,iLon] )

                # coef_center = -2 -(cosLatMHalf + cosLatPHalf) *cosLat 
                # coef_south  =   cosLatMHalf *cosLat
                # coef_north  =   cosLatPHalf *cosLat
                # coef_west   =  1 
                # coef_east   =  1
        # ## !!SOUTH POLE
        # iLat = 0
        # for iLon in range(self.nPhi):
        #     idx = iLat * self.nPhi + iLon
        #     coef_south_pole = - 4.0/ (self.R**2 * self.dLat**2)

        #     # screened poisson:   - (4/ (self.R**2 * (self.dLat/2)**2) + (2.0*OMEGA*self.sinLat0)**2/ (GRAVITY * self.d[iLat, iLon]))
        #     # scaled pole equation: -4*self.nPhi +(self.R**2 * self.dLat**2 *self.nPhi)* ((2.0*OMEGA*self.sinLat0)**2/(GRAVITY * self.d[iLat, iLon]))
            
        #     coef_ring = 4.0 / ( self.R**2 * self.dLat**2 * self.nPhi )
        #     # scaled pole equation: 4.0
            
        #     A[idx, idx] = coef_south_pole
        #     for jLon in range(self.nPhi):
        #         interior_idx = self.nPhi + jLon
        #         A[idx, interior_idx] = coef_ring

        # #!! NORTH POLE
        # iLat = self.nLat-1
        # for iLon in range(self.nPhi):
        #     idx = iLat * self.nPhi + iLon
        #     coef_north_pole = - 4.0/ (self.R**2 * self.dLat**2)
            
        #     # full screened poisson: -(4/(self.R**2 * (self.dLat/2)**2)+ (2.0*OMEGA*self.sinLat0)**2/(GRAVITY * self.d[iLat, iLon]))
        #     # scaled: -4*self.nPhi +(self.R**2 * (self.dLat/2)**2 *self.nPhi)* ((2.0*OMEGA*self.sinLat0)**2/(GRAVITY * self.d[iLat, iLon]))
            
        #     coef_ring = 4.0/(self.R**2 * self.dLat**2 * self.nPhi)
        #     # 4.0
        #     #
        #     A[idx, idx] = coef_north_pole
        #     # print('idx,idx: ', idx,idx)
        #     for jLon in range(self.nPhi):
        #         interior_idx = (self.nLat-2)*self.nPhi + jLon
        #         A[idx, interior_idx] = coef_ring
        #         # print('idx,interior_idx: ', idx,interior_idx)
        
        return A.tocsr()
    
    def Psi_build_rhs(self):
        print('APPENDING RHS VALUES FOR ψ')
        n_elements = self.nLat * self.nPhi
        RHS = np.zeros(n_elements)

        for iLat in range(self.nLat):
            Lat = self.Lat[iLat]
            cosLat = np.cos(Lat)
            sinLat = np.sin(Lat)
            for iLon in range(self.nPhi):
                idx = iLat * self.nPhi + iLon
                
                # COMMENT: SKIP LAND CELL and continue the nested loop
                if self.d[iLat,iLon] <= 0:
                    continue

                RHS[idx] = -  2.0 * self.OMEGA * sinLat * (R**2*self.dLat**2 *cosLat**2) * self.d[iLat,iLon]

        # #South pole (θ = -π/2 → sinθ = -1)
        # iLat = 0
        # for iLon in range(self.nPhi):
        #     idx = iLat * self.nPhi + iLon
        #     RHS[idx] = -2.0 * self.OMEGA * (-1.0)

        # # North pole (θ = +π/2 → sinθ = +1)
        # iLat = self.nLat - 1
        # for iLon in range(self.nPhi):
        #     idx = iLat * self.nPhi + iLon
        #     RHS[idx] = -2.0 * self.OMEGA * (1.0)

        return RHS
    
    @staticmethod
    def ADI_line_relaxation(A, b, nLat, nPhi, tol_rel=1e-8, maxiter=500, pre_smooth=True):
        N = A.shape[0]
        x = np.zeros(N)
        b_inf = np.linalg.norm(b, np.inf)
        if b_inf == 0:
            print('!!!!!   b_inf == 0')
            # b_inf = 1.0
        residuals = []

        A_csr = A.tocsr()
        diagA = A_csr.diagonal()
        M_inv = np.zeros_like(diagA)
        nonzero = diagA != 0
        M_inv[nonzero] = 1.0 / diagA[nonzero]

        for k in range(maxiter):
            # optional Jacobi pre-smoothing (acts like a simple preconditioner)
            if pre_smooth:
                r = b - A_csr.dot(x)
                x += M_inv * r

            # recompute residual
            r = b - A_csr.dot(x)

            #! Solving Rows: (each row iLat)
            for iLat in range(nLat):
                idx0 = iLat * nPhi
                idxs = np.arange(idx0, idx0 + nPhi)
                # build small block A_rr and rhs r_rr
                A_rr = A_csr[idxs, :][:, idxs].tocsc()
                r_rr = r[idxs]
                if np.linalg.norm(r_rr, np.inf) == 0:
                    continue
                # solve A_rr * delta = r_rr
                try:
                    delta = spsolve(A_rr, r_rr)
                except Exception:
                    # fallback to diagonal correction
                    delta = M_inv[idxs] * r_rr
                x[idxs] += delta

            # recompute residual after row sweeps
            r = b - A_csr.dot(x)

            #! Solving for columns: (each column iLon)
            for iLon in range(nPhi):
                idxs = np.arange(iLon, N, nPhi)  # pick every nPhi starting at iLon
                A_cc = A_csr[idxs, :][:, idxs].tocsc()
                r_cc = r[idxs]
                if np.linalg.norm(r_cc, np.inf) == 0:
                    continue
                try:
                    delta = spsolve(A_cc, r_cc)
                except Exception:
                    delta = M_inv[idxs] * r_cc
                x[idxs] += delta

            # compute relative inf-norm residual
            r = b - A_csr.dot(x)
            rel = np.linalg.norm(r, np.inf) / b_inf
            residuals.append(rel)
            
            if rel < tol_rel:
                break

        return x, residuals

    def PsiSolve(self):
        print('SOLVING SCREENED POISSON FOR ψ')
        A = self.Psi_build_matrix()
        
        print("Minimum diagonal: ", A.diagonal().min())
        print("MAX diagonal: ", A.diagonal().max())
        print("Norm: ", np.linalg.norm(A @ np.ones(A.shape[0])))

        RHS = self.Psi_build_rhs()

        #! 1. STABILIZATION - rescaling * or /
        # s = np.max(np.abs(A.diagonal()))
        # A   /=  s
        # RHS /= s

        # #! 2. Remove nullspace component
        # RHS -= np.mean(RHS)

        # #! 5. Optional tiny diagonal for CG
        # alpha = 1e-12
        # A.setdiag(A.diagonal() + alpha)

        N = A.shape[0]
        b_inf = np.linalg.norm(RHS, np.inf)
        if b_inf == 0:
            print('b_inf')
            b_inf = 1.0

        # Direct solver (reference)
        with tqdm(total=1, desc="Direct solver") as bar:
            PSI_direct = spsolve(A, RHS)
            bar.update(1)

        # Build ILU preconditioner
        print("Building ILU preconditioner...")
        ilu = spilu(
            A.tocsc(),
            permc_spec="COLAMD",
            diag_pivot_thresh=0.0,
            drop_tol=1e-4,
            fill_factor=10
        )
        M = LinearOperator(A.shape, matvec=lambda x: ilu.solve(x))

        maxiter = min(1000, N)
        tol_rel = 1e-8  

        # CG
        residuals_cg = []
        def cb_cg(xk):
            rinf = np.linalg.norm(RHS - A @ xk, np.inf) / b_inf
            residuals_cg.append(rinf)

        PSI_cg, info_cg = cg(A, RHS, M=M, rtol=tol_rel, maxiter=maxiter, callback=cb_cg)

        # BiCGSTAB preconditioned
        residuals_bicg = []
        def cb_bicg(xk):
            rinf = np.linalg.norm(RHS - A @ xk, np.inf) / b_inf
            residuals_bicg.append(rinf)

        PSI_bicg, info_bicg = bicgstab(A, RHS, M=M, rtol=tol_rel, maxiter=maxiter, callback=cb_bicg)
        print(f"CG info: {info_cg}")
        print(f"BiCGSTAB info: {info_bicg}")

        # --- ADI (line-relaxation) with Jacobi pre-smoother ---
        print('Running ADI (line-relaxation) solver...')
        PSI_adi, residuals_adi = self.ADI_line_relaxation(A=A, b=RHS, nLat=self.nLat, nPhi=self.nPhi, tol_rel=1e-8, maxiter=1000, pre_smooth=True)

        # --- Save residuals to CSV and plot ---
        def save_residuals(filename, res):
            with open(os.path.join(OUTPUT_PATH, filename), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['Iteration', 'RelativeInfResidual'])
                for i, r in enumerate(res, start=1):
                    writer.writerow([i, r])

        save_residuals('CG_residuals.csv', residuals_cg)
        save_residuals('BiCGSTAB_residuals.csv', residuals_bicg)
        save_residuals('ADI_residuals.csv', residuals_adi)

        # Plot residuals
        plt.figure()
        if residuals_cg:
            plt.semilogy(np.maximum(residuals_cg, 1e-16), label='CG (precond)')
        if residuals_bicg:
            plt.semilogy(np.maximum(residuals_bicg, 1e-16), label='BiCGSTAB (precond)')
        if residuals_adi:
            plt.semilogy(np.maximum(residuals_adi, 1e-16), label='ADI (line-relax)')
        plt.xlabel('Iteration')
        plt.ylabel('Relative inf-norm residual')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(OUTPUT_PATH, 'residuals_plot.png'), dpi=200)
        plt.close()

        print('Saved residuals and plot')

        # --- Save solutions (same as before) ---
        with open(os.path.join(OUTPUT_PATH,'b.csv'), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['iLat', 'iLon', 'b'])
            for iLat in range(self.nLat):
                for iLon in range(self.nPhi):
                    idx = iLat * self.nPhi + iLon
                    writer.writerow([iLat + 1, iLon + 1, RHS[idx]])

        # Save PSI arrays (reshape)
        PSI_cg = PSI_cg.reshape(self.nLat, self.nPhi)
        PSI_bicg = PSI_bicg.reshape(self.nLat, self.nPhi)
        PSI_adi = PSI_adi.reshape(self.nLat, self.nPhi)
        PSI_direct = PSI_direct.reshape(self.nLat, self.nPhi)

        # write CSVs for solutions
        def save_psi(filename, psi):
            with open(os.path.join(OUTPUT_PATH, filename), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['iLat', 'iLon', 'psi'])
                for iLat in range(self.nLat):
                    for iLon in range(self.nPhi):
                        idx = iLat * self.nPhi + iLon
                        writer.writerow([iLat + 1, iLon + 1, psi[iLat, iLon]])

        save_psi('PythonPsi_direct.csv', PSI_direct)
        save_psi('PythonPsi_CG.csv', PSI_cg)
        save_psi('PythonPsi_BICGSTAB.csv', PSI_bicg)
        save_psi('PythonPsi_ADI.csv', PSI_adi)

        # Save A triplets
        A_coo = A.tocoo()
        with open(os.path.join(OUTPUT_PATH,'AC.csv'), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['row', 'col', 'value'])
            for row, col, val in zip(A_coo.row + 1, A_coo.col + 1, A_coo.data):
                if val != 0:
                    writer.writerow([row, col, val])

        print('Saved AC.csv')

        return PSI_direct , PSI_cg, PSI_bicg , PSI_adi
    
    def Psi_compute_velocity(self, psi):
        print('COMPUTING THE VELOCITY FIELD FOR ψ')
        uθ = np.zeros_like(psi)     
        uλ = np.zeros_like(psi)

        #!COMMENT:  uθ = 1/(R cos θ)  ∂ψ/∂λ ;theta hat
        for iLat in range( self.nLat):               
            cosLat = np.cos(self.Lat[iLat])
            for iLon in range(self.nPhi):
                iLonM1 = iLon-1
                iLonP1 = iLon+1 
                if iLon == self.nPhi-1:
                    iLonP1 = 0
                if iLon == 0:
                    iLonM1 = self.nPhi-1

                uθ[iLat, iLon] = (psi[iLat, iLonP1] - psi[iLat, iLonM1]) / (2.0 * self.R * self.dLat * cosLat)
        #!COMMENT: BC ψ[i+1] - ψ[i] ≈ 0 which are already zeros; since zero array
        #!COMMENT:  uλ = -1/R  ∂ψ/∂θ ; lambda hat
        for iLat in range(self.nLat ):
            iLatP1 =  iLat+1 
            iLatM1 =  iLat-1

            if iLat == 0:
                iLatM1 = 0 
            if iLat == self.nLat-1:
                iLatP1 = self.nLat -1     

            for iLon in range(self.nPhi):
                uλ[iLat, iLon] = -(psi[iLatP1, iLon] - psi[iLatM1, iLon]) / (2.0 * self.R * self.dLat)
        
        #!COMMENT: BC ψ[i+1] - ψ[i] ≈ 0 at poles → uλ = 0 

        speed = np.sqrt(uθ **2 + uλ**2)

        print(f"ψ: uθ range:  [{uθ .min(): .4e}, {uθ .max(): .4e}]")
        print(f"ψ: uλ range: [{uλ.min(): .4e}, {uλ.max(): .4e}]")
        print(f"ψ: speed magnitude range: [{speed.min(): .4e}, {speed.max(): .4e}]")

        return uθ, uλ, speed
        
    def interactive_coefficient_matrix_plot(self,A,show_values=True):

        # dense matrix
        A_dense = A.toarray()
        nrows, ncols = A_dense.shape

        # mask zeros for transparency
        z = np.where(A_dense != 0, A_dense, np.nan)

        # hover text
        hover_text = [
            [f"({i},{j}): {A_dense[i,j]:.3e}" for j in range(ncols)]
            for i in range(nrows)
        ]

        fig = go.Figure(
            go.Heatmap(
                z=z,
                x=np.arange(ncols),
                y=np.arange(nrows),
                text=hover_text,
                hoverinfo='text',
                colorscale='Teal',
                zmin=np.min(A_dense),  # full matrix min
                zmax=np.max(A_dense),  # full matrix max
                colorbar=dict(
                    title="Value",
                    tickformat=".2e",     # scientific notation
                    ticks="outside",
                ),
            )
        )
        # add text annotations (optional)
        if show_values:
            for i in range(nrows):
                for j in range(ncols):
                    val = A_dense[i, j]
                    if val != 0:
                        fig.add_annotation(
                            x=j, y=i,
                            text=f"{val:.2g}",
                            showarrow=False,
                            font=dict(color='white', size=10),
                            xanchor='center', yanchor='middle'
                        )

        # layout with black theme
        fig.update_layout(
            template="plotly_dark",
            title=f"Coefficient Matrix (nLat={self.nLat}, nPhi={self.nPhi})",
            plot_bgcolor='black',
            paper_bgcolor='black',
            font_color='white',
            xaxis=dict(
                title="Column Index",
                tickmode='linear',
                dtick=1,
                side="top",
                showgrid=True,
                gridcolor="gray",
            ),
            yaxis=dict(
                title="Row Index",
                autorange='reversed',  # top = 0
                tickmode='linear',
                dtick=1,
                showgrid=True,
                gridcolor="gray",
            ),
            width=900 + ncols * 20,
            height=900 + nrows * 20,
        )

        fig.show()

    def plot_rhs_1d_plotly(self,RHS):
        idx = np.arange(len(RHS))

        fig = go.Figure(go.Scatter(
            x=idx,
            y=RHS,
            mode="markers+lines",
            marker=dict(size=6),
            hovertemplate="Index=%{x}<br>RHS=%{y:.3e}<extra></extra>"
        ))

        fig.update_layout(
            template="plotly_dark",
            title="RHS Vector (1D)",
            xaxis=dict(
                title="Linear index (0 → nLat·nPhi − 1)",
                tickmode="linear"
            ),
            yaxis=dict(
                title="RHS value",
                tickformat=".2e"
            ),
            plot_bgcolor="black",
            paper_bgcolor="black",
            font_color="white",
            width=900,
            height=400
        )

        fig.show()

    def plot_A_1d_plotly(self, A):
    
        A_coo = A.tocoo()

        idx = np.arange(len(A_coo.data))

        fig = go.Figure(go.Scatter(
            x=idx,
            y=A_coo.data,
            mode="markers",
            marker=dict(size=5),
            hovertemplate=(
                "nz-index=%{x}<br>"
                "A[%{customdata[0]}, %{customdata[1]}]=%{y:.3e}"
                "<extra></extra>"
            ),
            customdata=np.column_stack((A_coo.row, A_coo.col))
        ))

        fig.update_layout(
            template="plotly_dark",
            title="Coefficient Matrix A (1D view of nonzero entries)",
            xaxis=dict(title="Nonzero entry index"),
            yaxis=dict(title="Value", tickformat=".2e"),
            plot_bgcolor="black",
            paper_bgcolor="black",
            font_color="white",
            width=900,
            height=400
        )

        fig.show()

