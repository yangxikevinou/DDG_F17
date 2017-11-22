"use strict";

class Geometry {
	/**
	 * This class represents the geometry of a {@link module:Core.Mesh Mesh}. This includes information such
	 * as the position of vertices as well as methods to compute edge lengths, corner
	 * angles, face area, normals, discrete curvatures etc.
	 * @constructor module:Core.Geometry
	 * @param {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @param {module:LinearAlgebra.Vector[]} positions An array containing the position of each vertex in a mesh.
	 * @param {boolean} normalizePositions flag to indicate whether positions should be normalized. Default value is true.
	 * @property {module:Core.Mesh} mesh The mesh this class describes the geometry of.
	 * @property {module:LinearAlgebra.Vector[]} positions A dictionary mapping each vertex to a normalized position.
	 */
	constructor(mesh, positions, normalizePositions = true) {
		this.mesh = mesh;
		this.positions = {};
		for (let i = 0; i < positions.length; i++) {
			let v = this.mesh.vertices[i];
			let p = positions[i];

			this.positions[v] = p;
		}

		if (normalizePositions) {
			normalize(this.positions, mesh.vertices);
		}
	}

	/**
	 * Computes the vector along a halfedge.
	 * @method module:Core.Geometry#vector
	 * @param {module:Core.Halfedge} h The halfedge along which the vector needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vector(h) {
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];

		return b.minus(a);
	}

	/**
	 * Computes the length of an edge.
	 * @method module:Core.Geometry#length
	 * @param {module:Core.Edge} e The edge whose length needs to be computed.
	 * @returns {number}
	 */
	length(e) {
		return this.vector(e.halfedge).norm();
	}

	/**
	 * Computes the midpoint of an edge.
	 * @method module:Core.Geometry#midpoint
	 * @param {module:Core.Edge} e The edge whose midpoint needs to be computed.
	 * @returns {number}
	 */
	midpoint(e) {
		let h = e.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.twin.vertex];

		return (a.plus(b)).over(2);
	}

	/**
	 * Computes the mean edge length of all the edges in a mesh.
	 * @method module:Core.Geometry#meanEdgeLength
	 * @returns {number}
	 */
	meanEdgeLength() {
		let sum = 0;
		let edges = this.mesh.edges;
		for (let e of edges) {
			sum += this.length(e);
		}

		return sum / edges.length;
	}

	/**
	 * Computes the area of a face.
	 * @method module:Core.Geometry#area
	 * @param {module:Core.Face} f The face whose area needs to be computed.
	 * @returns {number}
	 */
	area(f) {
		if (f.isBoundaryLoop()) return 0.0;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return 0.5 * u.cross(v).norm();
	}

	/**
	 * Computes the total surface area of a mesh.
	 * @method module:Core.Geometry#totalArea
	 * @returns {number}
	 */
	totalArea() {
		let sum = 0.0;
		for (let f of this.mesh.faces) {
			sum += this.area(f);
		}

		return sum;
	}

	/**
	 * Computes the normal of a face.
	 * @method module:Core.Geometry#faceNormal
	 * @param {module:Core.Face} f The face whose normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	faceNormal(f) {
		if (f.isBoundaryLoop()) return undefined;

		let u = this.vector(f.halfedge);
		let v = this.vector(f.halfedge.prev).negated();

		return u.cross(v).unit();
	}

	/**
	 * Computes the centroid of a face.
	 * @method module:Core.Geometry#centroid
	 * @param {module:Core.Face} f The face whose centroid needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	centroid(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		return a.plus(b).plus(c).over(3);
	}

	/**
	 * Computes the circumcenter of a face.
	 * @method module:Core.Geometry#circumcenter
	 * @param {module:Core.Face} f The face whose circumcenter needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	circumcenter(f) {
		let h = f.halfedge;
		let a = this.positions[h.vertex];
		let b = this.positions[h.next.vertex];
		let c = this.positions[h.prev.vertex];

		if (f.isBoundaryLoop()) return a.plus(b).over(2);

		let ac = c.minus(a);
		let ab = b.minus(a);
		let w = ab.cross(ac);

		let u = (w.cross(ab)).times(ac.norm2());
		let v = (ac.cross(w)).times(ab.norm2());
		let x = (u.plus(v)).over(2 * w.norm2());

		return x.plus(a);
	}

	/**
	 * Computes an orthonormal bases for a face.
	 * @method module:Core.Geometry#orthonormalBases
	 * @param {module:Core.Face} f The face on which the orthonormal bases needs to be computed.
	 * @returns {module:LinearAlgebra.Vector[]} An array containing two orthonormal vectors tangent to the face.
	 */
	orthonormalBases(f) {
		let e1 = this.vector(f.halfedge).unit();

		let normal = this.faceNormal(f);
		let e2 = normal.cross(e1);

		return [e1, e2];
	}

	/**
	 * Computes the angle (in radians) at a corner.
	 * @method module:Core.Geometry#angle
	 * @param {module:Core.Corner} c The corner at which the angle needs to be computed.
	 * @returns {number} The angle clamped between 0 and π.
	 */
	angle(c) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		return Math.PI/2-Math.atan(this.cotan(c.halfedge));
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let u = this.vector(c.halfedge.prev).unit();
		let v = this.vector(c.halfedge.next).negated().unit();

		return Math.acos(Math.max(-1.0, Math.min(1.0, u.dot(v))));
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		return Math.PI/2-Math.atan(this.cotan(c.halfedge));
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Computes the cotangent of the angle opposite to a halfedge.
	 * @method module:Core.Geometry#cotan
	 * @param {module:Core.Halfedge} h The halfedge opposite to the angle whose cotangent needs to be computed.
	 * @returns {number}
	 */
	cotan(h) {
<<<<<<< HEAD
<<<<<<< HEAD
		let u=this.vector(h.next).negated();
		let v=this.vector(h.prev);
		if (h.next.next=h.prev) {
			return (u.dot(v))/((u.cross(v)).norm());
		}
		else {
			return undefined;
		}
<<<<<<< HEAD
=======
		if (h.onBoundary) return 0.0;

		let u = this.vector(h.prev);
		let v = this.vector(h.next).negated();

		return u.dot(v) / u.cross(v).norm();
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		if (h.onBoundary) return 0.0;

		let u = this.vector(h.prev);
		let v = this.vector(h.next).negated();

		return u.dot(v) / u.cross(v).norm();
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Computes the signed angle (in radians) between two adjacent faces.
	 * @method module:Core.Geometry#dihedralAngle
	 * @param {module:Core.Halfedge} h The halfedge (shared by the two adjacent faces) on which
	 * the dihedral angle is computed.
	 * @returns {number} The dihedral angle.
	 */
	dihedralAngle(h) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		let n1=this.faceNormal(h.face);
		let n2=this.faceNormal(h.twin.face);		
		return Math.atan2(this.vector(h).unit().dot(n1.cross(n2)),n1.dot(n2));
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		if (h.onBoundary || h.twin.onBoundary) return 0.0;

		let n1 = this.faceNormal(h.face);
		let n2 = this.faceNormal(h.twin.face);
		let w = this.vector(h).unit();

		let cosTheta = n1.dot(n2);
		let sinTheta = n1.cross(n2).dot(w);

		return Math.atan2(sinTheta, cosTheta);
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		let n1=this.faceNormal(h.face);
		let n2=this.faceNormal(h.twin.face);		
		return Math.atan2(this.vector(h).unit().dot(n1.cross(n2)),n1.dot(n2));
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Computes the barycentric dual area of a vertex.
	 * @method module:Core.Geometry#barycentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose barycentric dual area needs to be computed.
	 * @returns {number}
	 */
	barycentricDualArea(v) {
<<<<<<< HEAD
<<<<<<< HEAD
		let a=0.0;
		for (let f of v.adjacentFaces()){
			a+=this.area(f);
		}
		return a/3;
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let area = 0.0;
		for (let f of v.adjacentFaces()) {
			area += this.area(f) / 3;
		}

		return area;
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
	}

	/**
	 * Computes the circumcentric dual area of a vertex.
	 * @see {@link http://www.cs.cmu.edu/~kmcrane/Projects/Other/TriangleAreasCheatSheet.pdf}
	 * @method module:Core.Geometry#circumcentricDualArea
	 * @param {module:Core.Vertex} v The vertex whose circumcentric dual area needs to be computed.
	 * @returns {number}
	 */
	circumcentricDualArea(v) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let a=0.;
		for (let h of v.adjacentHalfedges()){
			a+=this.cotan(h)*Math.pow(this.length(h.edge),2);
			h=h.twin;
			a+=this.cotan(h)*Math.pow(this.length(h.edge),2);
		}
		return a/8;
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let area = 0.0;
		for (let h of v.adjacentHalfedges()) {
			let u2 = this.vector(h.prev).norm2();
			let v2 = this.vector(h).norm2();
			let cotAlpha = this.cotan(h.prev);
			let cotBeta = this.cotan(h);

			area += (u2 * cotAlpha + v2 * cotBeta) / 8;
		}

		return area;
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
	}

	/**
	 * Computes the normal at a vertex using the "equally weighted" method.
	 * @method module:Core.Geometry#vertexNormalEquallyWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalEquallyWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
			let normal = this.faceNormal(f);

			n.incrementBy(normal);
		}

		n.normalize();

		return n;
	}

	/**
	 * Computes the normal at a vertex using the "face area weights" method.
	 * @method module:Core.Geometry#vertexNormalAreaWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAreaWeighted(v) {
		let n = new Vector();
		for (let f of v.adjacentFaces()) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			n.incrementBy(this.faceNormal(f).times(this.area(f)));
		}
		n.normalize();
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
			let normal = this.faceNormal(f);
			let area = this.area(f);

			n.incrementBy(normal.times(area));
		}
<<<<<<< HEAD

		n.normalize();

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======

		n.normalize();

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
			n.incrementBy(this.faceNormal(f).times(this.area(f)));
		}
		n.normalize();
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		return n;
	}

	/**
	 * Computes the normal at a vertex using the "tip angle weights" method.
	 * @method module:Core.Geometry#vertexNormalAngleWeighted
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalAngleWeighted(v) {
		let n = new Vector();
		for (let c of v.adjacentCorners()) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			n.incrementBy(this.faceNormal(c.face).times(this.angle(c)));
		}
		n.normalize();
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
			let normal = this.faceNormal(c.halfedge.face);
			let angle = this.angle(c);

			n.incrementBy(normal.times(angle));
		}

		n.normalize();

<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
			n.incrementBy(this.faceNormal(c.face).times(this.angle(c)));
		}
		n.normalize();
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		return n;
	}

	/**
	 * Computes the normal at a vertex using the "gauss curvature" method.
	 * @method module:Core.Geometry#vertexNormalGaussCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalGaussCurvature(v) {
		let n = new Vector();
		for (let h of v.adjacentHalfedges()) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			if (!h.onBoundary) n.incrementBy(this.vector(h).unit().times(this.dihedralAngle(h)));
		}
		n.normalize();
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
			let weight = 0.5 * this.dihedralAngle(h) / this.length(h.edge);

			n.decrementBy(this.vector(h).times(weight));
		}

		n.normalize();

<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
			if (!h.onBoundary) n.incrementBy(this.vector(h).unit().times(this.dihedralAngle(h)));
		}
		n.normalize();
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		return n;
	}

	/**
	 * Computes the normal at a vertex using the "mean curvature" method (same as the "area gradient" method).
	 * @method module:Core.Geometry#vertexNormalMeanCurvature
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalMeanCurvature(v) {
		let n = new Vector();
		for (let h of v.adjacentHalfedges()) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			if (!h.onBoundary) n.incrementBy(this.vector(h).times(this.cotan(h)+this.cotan(h.twin)));
		}
		n.normalize();
=======
			let weight = 0.5 * (this.cotan(h) + this.cotan(h.twin));

			n.decrementBy(this.vector(h).times(weight));
		}

		n.normalize();

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
			let weight = 0.5 * (this.cotan(h) + this.cotan(h.twin));

			n.decrementBy(this.vector(h).times(weight));
		}

		n.normalize();

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
			if (!h.onBoundary) n.incrementBy(this.vector(h).times(this.cotan(h)+this.cotan(h.twin)));
		}
		n.normalize();
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		return n;
	}

	/**
	 * Computes the normal at a vertex using the "inscribed sphere" method.
	 * @method module:Core.Geometry#vertexNormalSphereInscribed
	 * @param {module:Core.Vertex} v The vertex on which the normal needs to be computed.
	 * @returns {module:LinearAlgebra.Vector}
	 */
	vertexNormalSphereInscribed(v) {
		let n = new Vector();
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		for (let h1 of v.adjacentHalfedges()) {
			let h2=h1.prev;
			if (!(h1.onBoundary && h2.onBoundary)) {
				h1=this.vector(h1);
				h2=this.vector(h2).negated();
				n.incrementBy(h1.cross(h2).over(Math.pow(h1.norm()*h2.norm(),2)));
			}
		}
		n.normalize();
<<<<<<< HEAD
=======
		for (let c of v.adjacentCorners()) {
			let u = this.vector(c.halfedge.prev);
			let v = this.vector(c.halfedge.next).negated();

			n.incrementBy(u.cross(v).over(u.norm2() * v.norm2()));
		}

		n.normalize();

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		for (let c of v.adjacentCorners()) {
			let u = this.vector(c.halfedge.prev);
			let v = this.vector(c.halfedge.next).negated();

			n.incrementBy(u.cross(v).over(u.norm2() * v.norm2()));
		}

		n.normalize();

>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		return n;
	}

	/**
	 * Computes the angle defect at a vertex (= 2π minus the sum of incident angles
	 * at an interior vertex or π minus the sum of incident angles at a boundary vertex).
	 * @method module:Core.Geometry#angleDefect
	 * @param {module:Core.Vertex} v The vertex whose angle defect needs to be computed.
	 * @returns {number}
	 */
	angleDefect(v) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let a=0.;
		for (let c of v.adjacentCorners()) a+=this.angle(c);
		if (v.onBoundary()) return Math.PI-a;
		else return 2.*Math.PI-a;
<<<<<<< HEAD
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let angleSum = 0.0;
		for (let c of v.adjacentCorners()) {
			angleSum += this.angle(c);
		}

		return v.onBoundary() ? Math.PI - angleSum : 2 * Math.PI - angleSum;
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Computes the (integrated) scalar gauss curvature at a vertex.
	 * @method module:Core.Geometry#scalarGaussCurvature
	 * @param {module:Core.Vertex} v The vertex whose gauss curvature needs to be computed.
	 * @returns {number}
	 */
	scalarGaussCurvature(v) {
		return this.angleDefect(v);
	}

	/**
	 * Computes the (integrated) scalar mean curvature at a vertex.
	 * @method module:Core.Geometry#scalarMeanCurvature
	 * @param {module:Core.Vertex} v The vertex whose mean curvature needs to be computed.
	 * @returns {number}
	 */
	scalarMeanCurvature(v) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		let H=0;
		for (let h of v.adjacentHalfedges()) H+=this.dihedralAngle(h)*this.length(h.edge);
		return H/2.;
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let sum = 0.0;
		for (let h of v.adjacentHalfedges()) {
			sum += 0.5 * this.length(h.edge) * this.dihedralAngle(h);
		}

		return sum;
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		let H=0;
		for (let h of v.adjacentHalfedges()) H+=this.dihedralAngle(h)*this.length(h.edge);
		return H/2.;
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Computes the total angle defect (= 2π times the euler characteristic of the mesh).
	 * @method module:Core.Geometry#totalAngleDefect
	 * @returns {number}
	 */
	totalAngleDefect() {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		return 2.*Math.PI*this.mesh.eulerCharacteristic();
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let totalDefect = 0.0;
		for (let v of this.mesh.vertices) {
			totalDefect += this.angleDefect(v);
		}

		return totalDefect;
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
		return 2.*Math.PI*this.mesh.eulerCharacteristic();
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
	 * @method module:Core.Geometry#principalCurvatures
	 * @param {module:Core.Vertex} v The vertex on which the principal curvatures need to be computed.
	 * @returns {number[]} An array containing the minimum and maximum principal curvature values at a vertex.
	 */
	principalCurvatures(v) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let H=this.scalarMeanCurvature(v)/this.circumcentricDualArea(v);
		let K=this.scalarGaussCurvature(v)/this.circumcentricDualArea(v);
		K=Math.sqrt(H*H-K);
		return [H-K, H+K];
<<<<<<< HEAD
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let A = this.circumcentricDualArea(v);
		let H = this.scalarMeanCurvature(v) / A;
		let K = this.angleDefect(v) / A;

		let discriminant = H * H - K;
		if (discriminant > 0) discriminant = Math.sqrt(discriminant);
		else discriminant = 0;
<<<<<<< HEAD

		let k1 = H - discriminant;
		let k2 = H + discriminant;

		return [k1, k2];
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======

		let k1 = H - discriminant;
		let k2 = H + discriminant;

		return [k1, k2];
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Builds a sparse laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#laplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	laplaceMatrix(vertexIndex) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let VS=this.mesh.vertices;
		let T=new Triplet(VS.length,VS.length);
		for (let v of VS) {
			let a=0.;
			for (let h of v.adjacentHalfedges()) {
				let b=0.5*(this.cotan(h)+this.cotan(h.twin));
				T.addEntry(-b,vertexIndex[v],vertexIndex[h.twin.vertex]);
				a+=b;
			}
			T.addEntry(a+1e-8,vertexIndex[v],vertexIndex[v]);
		}
<<<<<<< HEAD
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let V = this.mesh.vertices.length;
		let T = new Triplet(V, V);
		for (let v of this.mesh.vertices) {
			let i = vertexIndex[v];
			let sum = 1e-8;

			for (let h of v.adjacentHalfedges()) {
				let j = vertexIndex[h.twin.vertex];
				let weight = (this.cotan(h) + this.cotan(h.twin)) / 2;
				sum += weight;

				T.addEntry(-weight, i, j);
			}

			T.addEntry(sum, i, i);
		}

<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		return SparseMatrix.fromTriplet(T);
	}

	/**
	 * Builds a sparse diagonal mass matrix containing the barycentric dual area of each vertex
	 * of a mesh.
	 * @method module:Core.Geometry#massMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	massMatrix(vertexIndex) {
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
		let VS=this.mesh.vertices;
		let d=DenseMatrix.zeros(VS.length,1);
		for (let v of VS) d.set(this.barycentricDualArea(v),vertexIndex[v],0);
		return SparseMatrix.diag(d);
<<<<<<< HEAD
=======
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
		let V = this.mesh.vertices.length;
		let T = new Triplet(V, V);
		for (let v of this.mesh.vertices) {
			let i = vertexIndex[v];

			T.addEntry(this.barycentricDualArea(v), i, i);
		}

		return SparseMatrix.fromTriplet(T);
<<<<<<< HEAD
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 14617552d87fdb8f123aaad0ed286f6e1bc62ca5
=======
>>>>>>> 938cdf2889e021d8e4503139399a425641a2de14
	}

	/**
	 * Builds a sparse complex laplace matrix. The laplace operator is negative semidefinite;
	 * instead we build a positive definite matrix by multiplying the entries of the
	 * laplace matrix by -1 and shifting the diagonal elements by a small constant (e.g. 1e-8).
	 * @method module:Core.Geometry#complexLaplaceMatrix
	 * @param {Object} vertexIndex A dictionary mapping each vertex of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.ComplexSparseMatrix}
	 */
	complexLaplaceMatrix(vertexIndex) {
		let VS=this.mesh.vertices;
		let T=new ComplexTriplet(VS.length,VS.length);
		for (let v of VS) {
			let a=0.;
			for (let h of v.adjacentHalfedges()) {
				let b=0.5*(this.cotan(h)+this.cotan(h.twin));
				T.addEntry(new Complex(-b,0.),vertexIndex[v],vertexIndex[h.twin.vertex]);
				a+=b;
			}
			T.addEntry(new Complex(a+1e-8,0.),vertexIndex[v],vertexIndex[v]);
		}
		return ComplexSparseMatrix.fromTriplet(T);
	}
}

/**
 * Centers a mesh about the origin and rescales it to unit radius.
 * @global
 * @function module:Core.normalize
 * @param {module:LinearAlgebra.Vector[]} positions The position of each vertex in the vertices array.
 * @param {module:Core.Vertex[]} vertices The vertices of a mesh.
 * @param {boolean} rescale A flag indicating whether mesh positions should be scaled to a unit radius.
 */
function normalize(positions, vertices, rescale = true) {
	// compute center of mass
	let N = vertices.length;
	let cm = new Vector();
	for (let v of vertices) {
		let p = positions[v];

		cm.incrementBy(p);
	}
	cm.divideBy(N);

	// translate to origin and determine radius
	let radius = -1;
	for (let v of vertices) {
		let p = positions[v];

		p.decrementBy(cm);
		radius = Math.max(radius, p.norm());
	}

	// rescale to unit radius
	if (rescale) {
		for (let v of vertices) {
			let p = positions[v];

			p.divideBy(radius);
		}
	}
}
