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
        //注意这里面所有矩阵的命名规则是列行，比如A0叫VertexEdgeAdjacencyMatrix为Edge*Vertex

        /** Assigns indices to the input mesh's vertices, edges, and faces
         * @method module:Projects.SimplicialComplexOperators#assignElementIndices
         * @param {module:Core.Mesh} mesh The input mesh which we index.
         */
        assignElementIndices(mesh) {
                // TODO                
                mesh.indexElements();
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
                let vertexVector = DenseMatrix.zeros(this.mesh.vertices.length,1);
                //alert(vertexVector.get(0,0));
                let vertexIndex=indexElements(this.mesh.vertices);
                for (let vertices of subset.vertices)
                {
                        //vertexVector[this.mesh.vertices]
                        vertexVector.set(1, vertexIndex[vertices], 0);
                        //alert(vertexVector);
                }
                return vertexVector;
                
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
                let edgeVector=DenseMatrix.zeros(this.mesh.edges.length,1); 
                let edgeIndex=indexElements(this.mesh.edges);      
                for (let edges of subset.edges)
                {
                        edgeVector.set(1,edgeIndex[edges],0);
                        //alert(edgeVector.get(edges.index,0));
                }
                return edgeVector;
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
                let faceVector=DenseMatrix.zeros(this.mesh.faces.length,1);
                let faceIndex=indexElements(this.mesh.faces)
                for (let faces of subset.faces)
                {
                        faceVector.set(1,faceIndex[faces],0);
                }
                return faceVector;
        }
        /** Returns an array of the non-zero index of a vertor（this is a function I added)
        * @method module:Projects.SimplicialComplexOperators#nonZeroIndex
        * @param n X 1 denseMatrix
        * @returns an array containing 
        */
        nonZeroIndex(columnVector) {
                let nonZeroIndex=[];
                for (let idj=0;idj<columnVector.nRows();idj++) {
                        if (columnVector.get(idj,0)>0) {
                                nonZeroIndex.push(idj);
                                //alert(idj);                                
                        }
                }
                //alert(nonZeroIndex);
                return nonZeroIndex;

        }
        /** Returns the star of a subset.
         * @method module:Projects.SimplicialComplexOperators#star
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The star of the given subset.
         */
        star(subset) {
                // TODO
                //对于这个子集的所有点，线和面,找到包含它的所有单形
                let vertexVector=this.buildVertexVector(subset);
                let edgeVector=this.buildEdgeVector(subset);
                let faceVector=this.buildFaceVector(subset);
                // //jalert(vertexVector[vertexVector.length,1]);
                // let indexToVertex = new Map(); 
                // //let starset=MeshSubset.deepCopy(subset);
                let vertexAdjancentEdge=this.A0.timesDense(vertexVector);
                let vertexAdjancentFace=this.A1.timesDense(vertexAdjancentEdge);
                //let edgeAdjancentVertex=this.A0.transpose().timesDense(edgeVector);
                let edgeAdjancentFace=this.A1.timesDense(edgeVector);
                //let faceAdjancentVertex=(this.A1.timesSparse(this.A0)).transpose().timesDense(faceVector);
                //let faceAdjancentEdge=this.A1.transpose().timesDense(faceVector);
                //let starVerticesIndices=new Array();                
                //subset.addVertices(this.nonZeroIndex(edgeAdjancentVertex));
                //subset.addVertices(this.nonZeroIndex(faceAdjancentVertex));
                subset.addEdges(this.nonZeroIndex(vertexAdjancentEdge));
                //subset.addEdges(this.nonZeroIndex(faceAdjancentEdge));
                subset.addFaces(this.nonZeroIndex(vertexAdjancentFace));
                subset.addFaces(this.nonZeroIndex(edgeAdjancentFace));        

                return subset; // placeholder
        }

        /** Returns the closure of a subset.
         * @method module:Projects.SimplicialComplexOperators#closure
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The closure of the given subset.
         */
        closure(subset) {
                // TODO
                let vertexVector=this.buildVertexVector(subset);
                let edgeVector=this.buildEdgeVector(subset);
                let faceVector=this.buildFaceVector(subset);
                let vertexAdjancentEdge=this.A0.timesDense(vertexVector);
                let vertexAdjancentFace=this.A1.timesDense(vertexAdjancentEdge);
                let edgeAdjancentFace=this.A1.timesDense(edgeVector);

                let edgeAdjancentVertex=this.A0.transpose().timesDense(edgeVector);              

                let faceAdjancentVertex=(this.A1.timesSparse(this.A0)).transpose().timesDense(faceVector);
                let faceAdjancentEdge=this.A1.transpose().timesDense(faceVector);
                let addedFaceAdjancentVertex=(this.A1.timesSparse(this.A0)).transpose().timesDense(vertexAdjancentFace);
                let addedFaceAdjancentEdge=this.A1.transpose().timesDense(vertexAdjancentFace);

                //let vertexAdjacentVertex=edgeAdjancentVertex.timesDense(vertexAdjancentEdge);
                //let edgeAdjancentEdge=
                //let starVerticesIndices=new Array();                
                subset.addVertices(this.nonZeroIndex(edgeAdjancentVertex));
                subset.addVertices(this.nonZeroIndex(faceAdjancentVertex));
                //subset.addVertices(this.nonZeroIndex(addedFaceAdjancentVertex));

                //subset.addEdges(this.nonZeroIndex(vertexAdjancentEdge));
                subset.addEdges(this.nonZeroIndex(faceAdjancentEdge));
               // subset.addEdges(this.nonZeroIndex(addedFaceAdjancentEdge));
                //subset.addFaces(this.nonZeroIndex(vertexAdjancentFace));
                //subset.addFaces(this.nonZeroIndex(edgeAdjancentFace)); 
                return subset; // placeholder
        }

        /** Returns the link of a subset.
         * @method module:Projects.SimplicialComplexOperators#link
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {module:Core.MeshSubset} The link of the given subset.
         */
        link(subset) {
                // TODO
                let closureStarSet=MeshSubset.deepCopy(subset);
                let starClosureSet=MeshSubset.deepCopy(subset);                
                closureStarSet=this.closure(this.star(closureStarSet));
                starClosureSet=this.star(this.closure(starClosureSet));
                console.log(closureStarSet);
                let linkSet=MeshSubset.deepCopy(closureStarSet);
                linkSet.addSubset(starClosureSet);
                linkSet.deleteSubset(starClosureSet);
                //linkSet=linkSet.deleteSubset(starClosureSet);

//                 linkSet=linkSet.deleteSubset(starClosureSet);
                //subset=MeshSubset.deepCopy(linkSet);

                console.log(linkSet);
                return linkSet ; // placeholder
        }

        /** Returns true if the given subset is a subcomplex and false otherwise.
         * @method module:Projects.SimplicialComplexOperators#isComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {boolean} True if the given subset is a subcomplex and false otherwise.
         */
        isComplex(subset) {
                // TODO
                let closureSet=MeshSubset.deepCopy(subset);
                closureSet=this.closure(closureSet);
                if (subset.equals(closureSet)) {
                        return true;
                }
                else {return false;}
                
        }

        /** Returns the degree if the given subset is a pure subcomplex and -1 otherwise.
         * @method module:Projects.SimplicialComplexOperators#isPureComplex
         * @param {module:Core.MeshSubset} subset A subset of our mesh.
         * @returns {number} The degree of the given subset if it is a pure subcomplex and -1 otherwise.
         */
        isPureComplex(subset) {
                // TODO
                let vertexVector=this.buildVertexVector(subset);
                let edgeVector=this.buildEdgeVector(subset);
                let faceVector=this.buildFaceVector(subset);
                let faceAdjancentVertex=(this.A1.timesSparse(this.A0)).transpose().timesDense(faceVector);
                let edgeAdjancentVertex=this.A0.transpose().timesDense(edgeVector);                    
//                 let vertexMatrix = SparseMatrix.diag(vertexVector);
//                 let edgeMatrix=SparseMatrix.diag(edgeVector);
//                 let faceMatrix=SparseMatrix.diag(faceVector);
//                 let vertexEdgeMatrix=edgeMatrix.timesSparse(this.A0).timesSparse(vertexMatrix);
//                 let vertexFaceMatrix=faceMatrix.timesSparse(this.A1).timesSparse(this.A0).timesSparse(vertexMatrix);
                
//                 let vertexAdjancentEdge=this.A0.timesDense(vertexVector);
//                 let vertexAdjancentFace=this.A1.timesDense(vertexAdjancentEdge);
//                 let faceAdjancentVertex=(this.A1.timesSparse(this.A0)).transpose().timesDense(faceVector);
                if (this.isComplex(subset)) {
                        if (faceVector.sum()!=0) {
//                                 console.log(faceAdjancentVertex.sum());
//                                 console.log(vertexVector.sum());
                                for(let i=0;i<this.mesh.vertices.length;i++) {
                                        if (vertexVector.get(i,0)!=0) {
                                                if (faceAdjancentVertex.get(i,0)==0) {
                                                        return -1
                                                }        
                                        }
                                        if (faceAdjancentVertex.get(i,0)!=0) {
                                                if (vertexVector.get(i,0)==0) {
                                                        return -1
                                                }
                                        }
                                }
//                                 if (faceAdjancentVertex.sum()/2==vertexVector.sum()) {
//                                         return 2;        
                                                                
                                return 2;}
                                
                        //maybe below there is some problem
                        else if (edgeVector.sum()!=0) {
                                if (edgeAdjancentVertex.sum()==vertexVector.sum()) {
                                        return 1;
                                }
                                else {
                                        return -1;
                                }
                        }
                        else {return 0;}                                                                
                }
                else {return -1;}
        }

        /** Returns the boundary of a subset.
         * @method module:Projects.SimplicialComplexOperators#boundary
         * @param {module:Core.MeshSubset} subset A subset of our mesh. We assume <code>subset</code> is a pure subcomplex.
         * @returns {module:Core.MeshSubset} The boundary of the given pure subcomplex.
         */
        boundary(subset) {
                // TODO
                let boundarySet=new MeshSubset();
                let faceVector=this.buildFaceVector(subset);
                let edgeVector=this.buildEdgeVector(subset);
                //只有当只在faceEdgeNumVector中只出现一次的边为边界边
                if (this.isPureComplex(subset)==2) {
                        let faceEdgeNumVector=this.A1.transpose().timesDense(faceVector);
                        let boundaryEdgeVector=DenseMatrix.zeros(this.mesh.edges.length,1);
                        for(let i=0;i<this.mesh.edges.length;i++) {
                                if (faceEdgeNumVector.get(i,0)==1) {
                                        boundaryEdgeVector.set(1,i,0);
                                }
                        }
                        boundarySet.addEdges(this.nonZeroIndex(boundaryEdgeVector));
                         
                        return this.closure(boundarySet);
                }
                else if (this.isPureComplex(subset)==1) {
                        let edgeVertexNumVector=this.A0.transpose().timesDense(edgeVector);
                        let boundaryVertexVector=DenseMatrix.zeros(this.mesh.vertices.length,1);
                        for(let j=0;j<this.mesh.vertices.length;j++) {
                                if (edgeVertexNumVector.get(j,0)==1) {
                                        boundaryVertexVector.set(1,j,0);
                                }
                        }
                        boundarySet.addVertices(this.nonZeroIndex(boundaryVertexVector));
                        return boundarySet;
                }
                else if (this.isPureComplex(subset)==0) {
                        return boundarySet;
                }
                else {return boundarySet;}
                 // placeholder
        }
}
