"use strict";

class SpectralConformalParameterization {
	/**
	 * This class implements the {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf spectral conformal parameterization} algorithm to flatten
	 * surface meshes with boundaries conformally.
	 * @constructor module:Projects.SpectralConformalParameterization
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {Object} vertexIndex A dictionary mapping each vertex of the input mesh to a unique index.
	 */
	constructor(geometry) {
		this.geometry = geometry;
		this.vertexIndex = indexElements(geometry.mesh.vertices);
	}

	/**
	 * Builds the complex conformal energy matrix EC = ED - A.
	 * @private
	 * @method module:Projects.SpectralConformalParameterization#buildConformalEnergy
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	buildConformalEnergy() {
		let nV=this.geometry.mesh.vertices.length;
		let T=new ComplexTriplet(nV,nV);
		for (let f of this.geometry.mesh.boundaries) {
			for (let h of f.adjacentHalfedges()) {
				let i=this.vertexIndex[h.vertex];
				let j=this.vertexIndex[h.twin.vertex];
				T.addEntry(new Complex(0.,0.5),i,j);
				T.addEntry(new Complex(0.,-0.5),j,i);
			}
		}
		return this.geometry.complexLaplaceMatrix(this.vertexIndex).minus(ComplexSparseMatrix.fromTriplet(T)).timesComplex(new Complex(0.5,0.));
	}

	/**
	 * Flattens the input surface mesh with 1 or more boundaries conformally.
	 * @method module:Projects.SpectralConformalParameterization#flatten
	 * @returns {Object} A dictionary mapping each vertex to a vector of planar coordinates.
	 */
	flatten() {	
		let vertices = this.geometry.mesh.vertices;
		let flattening = this.geometry.positions;
		let xeig=Solvers.solveInversePowerMethod(this.buildConformalEnergy());

		for (let v of vertices) {
			let i = this.vertexIndex[v];
			let p = flattening[v];

			p.x = xeig.get(i, 0).re;
			p.y = xeig.get(i, 0).im;
			p.z = 0.;
		}

		// normalize flattening
		normalize(flattening, vertices);

		return flattening;
	}
}