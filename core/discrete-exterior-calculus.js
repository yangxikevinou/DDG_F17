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
<<<<<<< HEAD
<<<<<<< HEAD
		let VS=geometry.mesh.vertices;
		let d=DenseMatrix.zeros(VS.length,1);
		for (let v of VS){
			if (!v.onBoundary()) {
				d.set(geometry.barycentricDualArea(v),vertexIndex[v],0);
			}
		}
		return SparseMatrix.diag(d);
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let vertices = geometry.mesh.vertices;
		let V = vertices.length;
		let T = new Triplet(V, V);
		for (let v of vertices) {
			let i = vertexIndex[v];
			let area = geometry.barycentricDualArea(v);

			T.addEntry(area, i, i);
		}

		return SparseMatrix.fromTriplet(T);
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
<<<<<<< HEAD
<<<<<<< HEAD
		let ES=geometry.mesh.edges;
		let d=DenseMatrix.zeros(ES.length,1);
		for (let e of ES){
			if (!e.onBoundary()) {
				let h=e.halfedge;				
				d.set((geometry.cotan(h)+geometry.cotan(h.twin))/2,edgeIndex[e],0);
			}
		}
		return SparseMatrix.diag(d);
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let edges = geometry.mesh.edges;
		let E = edges.length;
		let T = new Triplet(E, E);
		for (let e of edges) {
			let i = edgeIndex[e];
			let w = (geometry.cotan(e.halfedge) + geometry.cotan(e.halfedge.twin)) / 2;

			T.addEntry(w, i, i);
		}

		return SparseMatrix.fromTriplet(T);
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
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
<<<<<<< HEAD
<<<<<<< HEAD
		let FS=geometry.mesh.faces;
		let d=DenseMatrix.zeros(FS.length,1);
		for (let f of FS){
			d.set(1/geometry.area(f),faceIndex[f],0);
		}
		return SparseMatrix.diag(d);
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let faces = geometry.mesh.faces;
		let F = faces.length;
		let T = new Triplet(F, F);
		for (let f of faces) {
			let i = faceIndex[f];
			let area = geometry.area(f);

			T.addEntry(1 / area, i, i);
		}

		return SparseMatrix.fromTriplet(T);
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
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
<<<<<<< HEAD
<<<<<<< HEAD
		let VS=geometry.mesh.vertices;
		let ES=geometry.mesh.edges;
		let T=new Triplet(ES.length,VS.length);
		for (let e of ES){
			let h=e.halfedge;
			T.addEntry(-1,edgeIndex[e],vertexIndex[h.vertex]);
			T.addEntry(1,edgeIndex[e],vertexIndex[h.twin.vertex]);
		}
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let edges = geometry.mesh.edges;
		let vertices = geometry.mesh.vertices;
		let E = edges.length;
		let V = vertices.length;
		let T = new Triplet(E, V);
		for (let e of edges) {
			let i = edgeIndex[e];
			let j = vertexIndex[e.halfedge.vertex];
			let k = vertexIndex[e.halfedge.twin.vertex];

			T.addEntry(1, i, j);
			T.addEntry(-1, i, k);
		}

<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
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
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let faces = geometry.mesh.faces;
		let edges = geometry.mesh.edges;
		let F = faces.length;
		let E = edges.length;
		let T = new Triplet(F, E);
		for (let f of faces) {
			let i = faceIndex[f];

			for (let h of f.adjacentHalfedges()) {
				let j = edgeIndex[h.edge];
				let sign = h.edge.halfedge === h ? 1 : -1;

				T.addEntry(sign, i, j);
			}
		}

<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		return SparseMatrix.fromTriplet(T);
	}
}
