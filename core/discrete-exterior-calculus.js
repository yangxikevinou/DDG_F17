"use strict";

/**
 * This class contains methods to build common {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
 * @memberof module:Core
 */
class DEC {
	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 0-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar0Form(geometry, vertexIndex) {
		let VS=geometry.mesh.vertices;
		let d=DenseMatrix.zeros(VS.length,1);
		for (let v of VS){
			if (!v.onBoundary()) {
				d.set(geometry.barycentricDualArea(v),vertexIndex[v],0);
			}
		}
		return SparseMatrix.diag(d);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		let ES=geometry.mesh.edges;
		let d=DenseMatrix.zeros(ES.length,1);
		for (let e of ES){
			if (!e.onBoundary()) {
				let h=e.halfedge;				
				d.set((geometry.cotan(h)+geometry.cotan(h.twin))/2,edgeIndex[e],0);
			}
		}
		return SparseMatrix.diag(d);
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 2-forms.
	 * By convention, the area of a vertex is 1.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar2Form(geometry, faceIndex) {
		let FS=geometry.mesh.faces;
		let d=DenseMatrix.zeros(FS.length,1);
		for (let f of FS){
			d.set(1/geometry.area(f),faceIndex[f],0);
		}
		return SparseMatrix.diag(d);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 0-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative0Form(geometry, edgeIndex, vertexIndex) {
		let VS=geometry.mesh.vertices;
		let ES=geometry.mesh.edges;
		let T=new Triplet(ES.length,VS.length);
		for (let e of ES){
			let h=e.halfedge;
			T.addEntry(-1,edgeIndex[e],vertexIndex[h.vertex]);
			T.addEntry(1,edgeIndex[e],vertexIndex[h.twin.vertex]);
		}
		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse matrix encoding the exterior derivative on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} faceIndex A dictionary mapping each face of a mesh to a unique index.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildExteriorDerivative1Form(geometry, faceIndex, edgeIndex) {
		let ES=geometry.mesh.edges;
		let FS=geometry.mesh.faces;
		let T=new Triplet(FS.length,ES.length);
		for (let f of FS){
			for (let h of f.adjacentHalfedges()){
				let e=h.edge;
				if (h==e.halfedge) T.addEntry(1,faceIndex[f],edgeIndex[e]);
				else T.addEntry(-1,faceIndex[f],edgeIndex[e]);
			}
		}
		return SparseMatrix.fromTriplet(T);
	}
}
