"use strict";

/**
 * This class contains methods to build common {@link https://cs.cmu.edu/~kmcrane/Projects/DDG/paper.pdf discrete exterior calculus} operators.
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
		// TODO
		let hodgeStar0Triplet=new Triplet(geometry.mesh.vertices.length,geometry.mesh.vertices.length);
		for (let v of geometry.mesh.vertices) {
			hodgeStar0Triplet.addEntry(geometry.barycentricDualArea(v),vertexIndex[v],vertexIndex[v]);
		}
		return SparseMatrix.fromTriplet(hodgeStar0Triplet); // placeholder
	}

	/**
	 * Builds a sparse diagonal matrix encoding the Hodge operator on 1-forms.
	 * @static
	 * @param {module:Core.Geometry} geometry The geometry of a mesh.
	 * @param {Object} edgeIndex A dictionary mapping each edge of a mesh to a unique index.
	 * @returns {module:LinearAlgebra.SparseMatrix}
	 */
	static buildHodgeStar1Form(geometry, edgeIndex) {
		// TODO
		let hodgeStar1Triplet=new Triplet(geometry.mesh.edges.length,geometry.mesh.edges.length);
		for (let e of geometry.mesh.edges) {			
			hodgeStar1Triplet.addEntry(0.5*(geometry.cotan(e.halfedge)+geometry.cotan(e.halfedge.twin)),edgeIndex[e],edgeIndex[e]);
		}
		return SparseMatrix.fromTriplet(hodgeStar1Triplet); // placeholder
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
		// TODO
		let hodgeStar2Triplet=new Triplet(geometry.mesh.faces.length,geometry.mesh.faces.length);
		for (let f of geometry.mesh.faces) {
			hodgeStar2Triplet.addEntry(1/geometry.area(f),faceIndex[f],faceIndex[f]);
		}

		return SparseMatrix.fromTriplet(hodgeStar2Triplet); // placeholder
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
		// TODO
		let exteriorDerivative0Triplet=new Triplet(geometry.mesh.edges.length,geometry.mesh.vertices.length);
		for (let e of geometry.mesh.edges) {
			exteriorDerivative0Triplet.addEntry(-1,edgeIndex[e],vertexIndex[e.halfedge.vertex]);
			exteriorDerivative0Triplet.addEntry(1,edgeIndex[e],vertexIndex[e.halfedge.next.vertex])
		} 
		return SparseMatrix.fromTriplet(exteriorDerivative0Triplet); // placeholder
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
		// TODO
		let exteriorDerivative1Triplet=new Triplet(geometry.mesh.faces.length,geometry.mesh.edges.length);
		for (let f of geometry.mesh.faces) {
			for (let h of f.adjacentHalfedges()) {
				if (h.edge.halfedge===h) {
					exteriorDerivative1Triplet.addEntry(1,faceIndex[f],edgeIndex[h.edge]);	
				}
				else if (h.edge.halfedge.twin===h) {
					exteriorDerivative1Triplet.addEntry(-1,faceIndex[f],edgeIndex[h.edge]);	
				}			
			}
		}
		return SparseMatrix.fromTriplet(exteriorDerivative1Triplet); // placeholder
	}
}
