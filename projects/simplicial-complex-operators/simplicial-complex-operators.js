"use strict";

/**
 * @module Projects
 */
class SimplicialComplexOperators {

        /** This class implements various operators (e.g. boundary, star, link) on a mesh.
         * @constructor module:Projects.SimplicialComplexOperators
         * @param {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:Core.Mesh} mesh The input mesh this class acts on.
         * @property {module:LinearAlgebra.SparseMatrix} A0 The vertex-edge adjacency matrix of <code>mesh</code>.
         * @property {module:LinearAlgebra.SparseMatrix} A1 The edge-face adjacency matrix of <code>mesh</code>.
         */
        constructor(mesh) {
                this.mesh = mesh;
                this.assignElementIndices(this.mesh);

                this.A0 = this.buildVertexEdgeAdjacencyMatrix(this.mesh);
                this.A1 = this.buildEdgeFaceAdjacencyMatrix(this.mesh);
        }

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                // TODO
                mesh.indexElements(mesh);
                // let vertexIndex={};
                // let edgeIndex={};
                // let faceIndex={};
                // let i=0;
                // let j=0;
                // let k=0;
                // this.mesh=mesh;
                // for(let vertices of this.mesh.vertices)
                // {
                //         vertexIndex[vertices]=i++;
                // }
                // for(let edges of this.mesh.edges)
                // {
                //         edgeIndex[edges]=j++;
                // }
                // for(let faces of this.mesh.faces)
                // {
                //         faceIndex[faces]=k++;
                // }
        }

        /** Returns the vertex-edge adjacency matrix of the given mesh.
         * @method module:Projects.SimplicialComplexOperators#buildVertexEdgeAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.DenseMatrix} The vertex-edge adjacency matrix of the given mesh.
         */
        buildVertexEdgeAdjacencyMatrix(mesh) {
                // TODO
                let edgeVertexTriplet = new Triplet(mesh.edges.length,mesh.vertices.length);
                for(let edge of this.mesh.edges)
                {
                        edgeVertexTriplet.addEntry(1,edge.index,edge.halfedge.vertex.index);
                        edgeVertexTriplet.addEntry(1,edge.index,edge.halfedge.twin.vertex.index);
                };
                let edgeVertexMatrix=SparseMatrix.fromTriplet(edgeVertexTriplet);
                return edgeVertexMatrix;
        }       

        /** Returns the edge-face adjacency matrix.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeFaceAdjacencyMatrix
         * @param {module:Core.Mesh} mesh The mesh whose adjacency matrix we compute.
         * @returns {module:LinearAlgebra.DenseMatrix} The edge-face adjacency matrix of the given mesh.
         */
        buildEdgeFaceAdjacencyMatrix(mesh) {
                // TODO
                let faceEdgeTriplet = new Triplet(mesh.faces.length,mesh.edges.length);
                for(let faces of this.mesh.faces)
                {
                        for (let edges of faces.adjacentEdges())
                        {
                                faceEdgeTriplet.addEntry(1,faces.index,edges.index,);                                
                        }
                        
                }
                let faceEdgeMatrix=SparseMatrix.fromTriplet(faceEdgeTriplet)
                return faceEdgeMatrix;
        }

        /** Returns a column vector representing the vertices of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildVertexVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |V| entries. The ith entry is 1 if
         *  vertex i is in the given subset and 0 otherwise
         */
        buildVertexVector(subset) {
                // TODO
                let V=DenseMatrix.zeros(mesh.vertices.length,1)
                for (let vertices of this.subset.vertices)
                {
                        V[vertices.index-1]=1;
                }
                return V;
                
        }

        /** Returns a column vector representing the edges of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildEdgeVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |E| entries. The ith entry is 1 if
         *  edge i is in the given subset and 0 otherwise
         */
        buildEdgeVector(subset) {
                // TODO
                let V=DenseMatrix.zeros(mesh.edges.length,1)
                for (let edges of this.subset.edges)
                {
                        V[edges.index-1]=1;
                }
                return V;
        }

        /** Returns a column vector representing the faces of the
         * given subset.
         * @method module:Projects.SimplicialComplexOperators#buildFaceVector
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:LinearAlgebra.DenseMatrix} A column vector with |F| entries. The ith entry is 1 if
         *  face i is in the given subset and 0 otherwise
         */
        buildFaceVector(subset) {
                // TODO
                let V=DenseMatrix.zeros(mesh.faces.length,1)
                for (let faces of this.subset.faces)
                {
                        V[faces.index-1]=1;
                }
                return V;
        }

        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO
                //对于这个子集的所有点，线和面,找到包含它的所有单形
                let starset=MeshSubset.deepCopy(subset);
                for (let vertex of subset.vertices)
                {
                        starset.addEdges(vertex.adjacentEdges)
                        starset.addFaces(vertex.adjacentFaces)
                        //find indexmatrix
                        // A0[vertices.index]
                        // starset.addEdges(edges[index]);
                } ;
                for (let face of subset.faces)
                {
                        //starset.addEdges(edge.vertex)
                        starset.addFaces(face.adjacentFaces)
                        //find indexmatrix
                        // A0[vertices.index]
                        // starset.addEdges(edges[index]);
                } ;
                return starset; // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO

                return subset; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO

                return subset; // placeholder
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // TODO
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // TODO
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // TODO

                return subset; // placeholder
        }
}
