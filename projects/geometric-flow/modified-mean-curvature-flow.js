"use strict";

class ModifiedMeanCurvatureFlow extends MeanCurvatureFlow {
	/**
	 * This class performs a {@link http://www.cs.jhu.edu/~misha/MyPapers/SGP12.pdf modified version} of {@link https://www.cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf mean curvature flow} on a surface mesh.
	 * @constructor module:Projects.ModifiedMeanCurvatureFlow
	 * @augments module:Projects.MeanCurvatureFlow
	 * @param {module:Core.Geometry} geometry The input geometry of the mesh this class acts on.
	 * @property {module:LinearAlgebra.SparseMatrix} A The laplace matrix of the input mesh.
	 */
	constructor(geometry) {
		super(geometry);
<<<<<<< HEAD
<<<<<<< HEAD

		// build the laplace matrix
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		this.A = geometry.laplaceMatrix(this.vertexIndex);
	}

	/**
	 * @inheritdoc
	 */
	buildFlowOperator(M, h) {
<<<<<<< HEAD
<<<<<<< HEAD
		return this.A.timesReal(h).plus(M);
=======
		// F = M + hA
		return M.plus(this.A.timesReal(h));
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		// F = M + hA
		return M.plus(this.A.timesReal(h));
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
	}
}
