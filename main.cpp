#include <igl/boundary_loop.h>
#include <igl/cat.h>
#include <igl/colon.h>
#include <igl/slice.h>
#include <igl/upsample.h>
#include <igl/false_barycentric_subdivision.h>
#include <igl/decimate.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/harmonic.h>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>

using Eigen::VectorXi;
using Eigen::MatrixXd;
using Eigen::MatrixXi;
using Eigen::VectorXd;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;

void printHelpExit() {
	printf("Invalid command line arguments specified!\n\n");

	printf("USAGE: hole_fixer [options]\n\n");

	printf("OPTIONS: \n");

	printf("  -in\t\t\ttarget mesh file in .off-format, with a hole\n");
	printf("  -out\t\t\toutput mesh file in .off-format\n");
	printf("  -outfaces\t\tHow many faces to decimate the mesh to\n");
	printf("  -upsample\t\tHow much upsampling to use when creating the patch\n");

	exit(1);
}

const char* findToken(const char* param, int argc, char* argv[]) {
	const char* token = nullptr;
	for (int i = 0; i < argc; ++i) {
		if (strcmp(argv[i], param) == 0) {
			if (i + 1 < argc) {
				token = argv[i + 1];
				break;
			}
		}
	}

	if (token == nullptr) {
		printf("Could not find command-line parameter %s\n", param);
		return nullptr;
	}

	return token;
}

const char* parseStringParam(const char* param, int argc, char* argv[]) {
	const char* token = findToken(param, argc, argv);
	return token;
}

bool parseIntParam(const char* param, int argc, char* argv[], unsigned int& out) {
	const char* token = findToken(param, argc, argv);
	if (token == nullptr)
		return false;

	int r = sscanf(token, "%u,", &out);
	if (r != 1 || r == EOF) {
		return false;
	}
	else {
		return true;
	}
}

void get2ringBoundaryMesh(const MatrixXi& originalF, MatrixXi& nRingF, VectorXi& originalLoop)
{
	int m = originalF.array().maxCoeff()+1;
	
	igl::boundary_loop(originalF, originalLoop);	

	Eigen::SparseVector<int> isBoundary(m);
	for (int i = 0; i < originalLoop.size(); i++)
	{
		isBoundary.coeffRef(originalLoop(i)) = 1;
	}

	Eigen::SparseMatrix<int> adjMatrix(m, m);
	igl::adjacency_matrix(originalF, adjMatrix);

	Eigen::SparseVector<int> n_ring_boundary = (adjMatrix)* isBoundary;

	//Triangles in n-ring
	//A face that belongs to 2-ring boundary contains at least 1 vertex from boundary loop
	int nonZeros = n_ring_boundary.nonZeros();
	nRingF = MatrixXi(nonZeros, 3);
	int j = 0;
	for (int i = 0; i < originalF.rows(); i++)
	{
		int a = originalF(i, 0);
		int b = originalF(i, 1);
		int c = originalF(i, 2);

		if (isBoundary.coeffRef(a) > 0 || isBoundary.coeffRef(b) > 0 || isBoundary.coeffRef(c) > 0)
		{
			nRingF(j, 0) = a;
			nRingF(j, 1) = b;
			nRingF(j, 2) = c;
			j++;
		}
	}
}

void subdivideInnerRingBoundary(
	const MatrixXd& originalV, 
	const MatrixXi& originalF, 
	const MatrixXi& ringF,
	const VectorXi& ringBoundary,
	int divisions,
	MatrixXd& ringDivV,
	MatrixXi& ringDivF
)
{
	if(divisions<=2)
	{
		ringDivF = ringF;
		// Has more vertices than we need 
		// ( should be a one liner somwhere to remove the not needed ones
		ringDivV = originalV; 
		return;
	}
	


	Eigen::SparseVector<int> isBoundary((int)originalV.rows());
	for (int i = 0; i < ringBoundary.size(); i++)
	{
		isBoundary.coeffRef(ringBoundary(i)) = 1;
	}

	std::cout << "Rows before" << originalV.rows();
	//upsampling 
	int nonZeros = isBoundary.nonZeros();
	int sizeIncr = nonZeros*(divisions-2);

	ringDivV = originalV;
	ringDivV.conservativeResize(ringDivV.rows()+sizeIncr, 3);

	int counterV = originalV.rows();

	

	ringDivF = ringF;
    int counterF = ringF.rows();
	ringDivF.conservativeResize(ringDivF.rows() + sizeIncr, 3);

	int initF = counterF;
	for (int i = 0; i < initF; i++)
	{
		int a = ringDivF(i, 0);
		int b = ringDivF(i, 1);
		int c = ringDivF(i, 2);
		
		auto boundaryCounter = [&](int i){
		    return (int)isBoundary.coeffRef(i)>0;
		};
		int boundaryCount = boundaryCounter(a) + boundaryCounter(b) + boundaryCounter(c);

		if (boundaryCount==2)
		{
						
			int v[2]={0,0},vi=0;
			if(boundaryCounter(a))
				v[vi++]=a;
			if(boundaryCounter(b))
				v[vi++]=b;
			if(boundaryCounter(c))
				v[vi++]=c;

			//auto divisions = 3;

			auto midPoint = [&](int i) -> VectorXd {
				double s0 = ((double)divisions - i - 1)/((double)divisions-1);
				double s1 = 1-s0;
				return s0 * ringDivV.row(v[0]) + s1*ringDivV.row(v[1]);	
			};

			for (int i=1;i<=divisions-2;i++)
			{
			    ringDivV.row(counterV+i-1) = midPoint(i);
			}
		 
			{
				std::vector<int> tri = {a,b,c};
				if(!isBoundary.coeffRef(b)){
					std::rotate(tri.begin(), tri.begin()+1,tri.end());
				}else if(!isBoundary.coeffRef(c)){
					std::rotate(tri.begin(), tri.begin()+2,tri.end());
				}
				a = tri[0];
				b = tri[1];
				c = tri[2];
			}			

			// Replace the first big triangle with the first subdivision
			ringDivF(i, 0) = a;
			ringDivF(i, 1) = b;
			ringDivF(i, 2) = counterV;

			for (int i = 1; i < divisions - 2; i++)
			{
				auto j = counterF + i - 1;
				ringDivF(j, 0) = a;
				ringDivF(j, 1) = counterV + i - 1;
				ringDivF(j, 2) = counterV + i;
			}

			{
				auto j = counterF + divisions - 3;
				ringDivF(j, 0) = a;
				ringDivF(j, 1) = counterV + divisions - 3;
				ringDivF(j, 2) = c;
			}
						
			
			counterF+=divisions-2;
			counterV+=divisions-2;


		}

	}
	std::cout << "Rows after" << originalV.rows();
}

void createMeshPatch(const VectorXi& originalLoop, const MatrixXd& originalV, MatrixXd& patchV, MatrixXi& patchF)
{
	// compute boundary center.
	VectorXd bcenter(3);
	{
		VectorXi R = originalLoop;
		VectorXi C(3); C << 0, 1, 2;
		MatrixXd B;
		MatrixXd A = originalV;
		igl::slice(A, R, C, B);
		bcenter = (1.0f / originalLoop.size()) * B.colwise().sum();
	}

	// a flat patch that fills the hole.
	//patchV = MatrixXd(originalLoop.size() + 1, 3); // patch will have an extra vertex for the center vertex.
	//patchF = MatrixXi(originalLoop.size(), 3);

	{
		VectorXi R = originalLoop;
		VectorXi C(3); C << 0, 1, 2;
		MatrixXd A = originalV;
		MatrixXd temp1;
		igl::slice(A, R, C, temp1);

		MatrixXd temp2(1, 3);
		temp2 << bcenter(0), bcenter(1), bcenter(2);

		// patch vertices will be the boundary vertices, plus a center vertex. concat these together.
		igl::cat(1, temp1, temp2, patchV);

		// create triangles that connect boundary vertices and the center vertex.
		for (int i = 0; i < originalLoop.size(); ++i) {
			patchF(i, 2) = (int)originalLoop.size();
			patchF(i, 1) = i;
			patchF(i, 0) = (1 + i) % originalLoop.size();
		}

		// also upsample patch. patch and original mesh will have the same upsampling level now
		// making it trivial to fuse them together.		

		//igl::upsample(Eigen::MatrixXd(patchV), Eigen::MatrixXi(patchF), patchV, patchF, upsampleN);
	}
}

void fuseMeshes(const MatrixXd& patchV, const MatrixXi& patchF, const MatrixXd& originalV, const MatrixXi& originalF, MatrixXd& fairedV, MatrixXi& fairedF)
{
	// the final mesh, where the patch and original mesh has been fused together.
	std::vector<std::vector<double>> fusedV;
	std::vector<std::vector<int>> fusedF;

	int index = 0; // vertex index counter for the fused mesh.

				   // first, we add the upsampled patch to the fused mesh.
	{
		for (; index < patchV.rows(); ++index) {
			fusedV.push_back({ patchV(index, 0), patchV(index, 1), patchV(index, 2) });
		}

		int findex = 0;
		for (; findex < patchF.rows(); ++findex) {
			fusedF.push_back({ patchF(findex, 0), patchF(findex, 1), patchF(findex, 2) });
		}
	}

	// now, we fuse the patch and the original mesh together.
	{
		// remaps indices from the original mesh to the fused mesh.
		std::map<int, int> originalToFusedMap;

		for (int itri = 0; itri < originalF.rows(); ++itri) {

			int triIndices[3];
			for (int iv = 0; iv < 3; ++iv) {

				int vertexIndex = originalF(itri, iv);

				int ret;

				if (originalToFusedMap.count(vertexIndex) == 0) {
					int foundMatch = -1;

					// the vertices at the boundary are the same, for both the two meshes(patch and original mesh).
					// we ensure that these vertices are not duplicated.
					// this is also how we ensure that the two meshes are fused together.
					for (int jj = 0; jj < patchV.rows(); ++jj) {
						VectorXd u(3); u << fusedV[jj][0], fusedV[jj][1], fusedV[jj][2];
						VectorXd v(3); v << originalV(vertexIndex, 0), originalV(vertexIndex, 1), originalV(vertexIndex, 2);

						if ((u - v).norm() < 0.00001) {
							foundMatch = jj;
							break;
						}
					}

					if (foundMatch != -1) {
						originalToFusedMap[vertexIndex] = foundMatch;
						ret = foundMatch;
					}
					else {
						fusedV.push_back({ originalV(vertexIndex, 0), originalV(vertexIndex, 1), originalV(vertexIndex, 2) });
						originalToFusedMap[vertexIndex] = index;
						ret = index;
						index++;
					}
				}
				else {
					ret = originalToFusedMap[vertexIndex];
				}

				triIndices[iv] = ret;
			}

			fusedF.push_back({
				triIndices[0],
				triIndices[1],
				triIndices[2] });


		}

	}

	fairedV = MatrixXd(fusedV.size(), 3); //vertices
	fairedF = MatrixXi (fusedF.size(), 3); //faces


	for (int vindex = 0; vindex < fusedV.size(); ++vindex) {
		auto r = fusedV[vindex];

		fairedV(vindex, 0) = r[0];
		fairedV(vindex, 1) = r[1];
		fairedV(vindex, 2) = r[2];
	}

	for (int findex = 0; findex < fusedF.size(); ++findex) {
		auto r = fusedF[findex];

		fairedF(findex, 0) = r[0];
		fairedF(findex, 1) = r[1];
		fairedF(findex, 2) = r[2];
	}

}

void  smoothMeshWithFixed2ring(MatrixXd& meshV, MatrixXi& meshF, MatrixXd& smoothedMeshV)

{	
	int k = 2;
	MatrixXi ringF; 
	VectorXi ringBoundary;

	get2ringBoundaryMesh(meshF, ringF, ringBoundary);
	// surface fairing simply means that we solve the equation
	// Delta^2 f = 0
	// with appropriate boundary conditions.
	// this function igl::harmonic from libigl takes care of that.
	std::set<int> bs; 
	for(int j = 0; j < ringF.array().size();j++){
	  bs.insert(ringF.array()(j));
	}
	MatrixXi b(bs.size(),1);
	int j = 0;
	for(int idx:bs)
	   b(j++)=idx;
	
	// Boundary condition values
	MatrixXd bc(b.rows(),3);
	
	for(int j = 0 ; j<b.rows();j++){
		bc.row(j) = meshV.row(b(j));
	}

	igl::harmonic(meshV, meshF, b, bc, k, smoothedMeshV);
}

void  remove2ringBoundary(const MatrixXi& originalF, MatrixXi& nRingF)
{
	int m = originalF.array().maxCoeff() + 1;
	VectorXi originalLoop; // indices of the boundary of the hole. 
	igl::boundary_loop(originalF, originalLoop);

	Eigen::SparseVector<int> isBoundary(m);
	for (int i = 0; i < originalLoop.size(); i++)
	{
		isBoundary.coeffRef(originalLoop(i)) = 1;
	}

	Eigen::SparseMatrix<int> adjMatrix(m, m);
	igl::adjacency_matrix(originalF, adjMatrix);

	Eigen::SparseVector<int> n_ring_boundary = (adjMatrix)*isBoundary;

	//Triangles in original mesh without 2-ring boundary
	//Add only those faces, none of the vertices of which belong to the boundary
	int nonZeros = n_ring_boundary.nonZeros();
	nRingF = MatrixXi(originalF.rows()-nonZeros, 3);
	int j = 0;
	for (int i = 0; i < originalF.rows(); i++)
	{
		int a = originalF(i, 0);
		int b = originalF(i, 1);
		int c = originalF(i, 2);

		if ((isBoundary.coeffRef(a) == 0 && isBoundary.coeffRef(b) == 0 && isBoundary.coeffRef(c) == 0))
		{
			nRingF(j, 0) = a;
			nRingF(j, 1) = b;
			nRingF(j, 2) = c;
			j++;
		}
	}
}

bool fillHole(const MatrixXd& originalV, const MatrixXi& originalF, int upsampleN, MatrixXd& fairedV, MatrixXi& fairedF)
{
	VectorXi originalLoop; // indices of the boundary of the hole. 
	igl::boundary_loop(originalF, originalLoop);

	if (originalLoop.size() == 0) {
		fairedV = originalV;
		fairedF = originalF;
		return false;
	}
	MatrixXi nRingF;
	VectorXi nRingBoundary;
	//Create N-ring boundary face matrix
	get2ringBoundaryMesh(originalF, nRingF, nRingBoundary);
	igl::writeOFF("2ring.off", originalV, nRingF);

	//Subdivide inner ring boundary
	MatrixXd nRingDivV;
	MatrixXi nRingDivF;
	int divisions = pow(2, upsampleN) + 1;
	subdivideInnerRingBoundary(originalV, originalF, nRingF, nRingBoundary, divisions, nRingDivV, nRingDivF);
	igl::writeOFF("2ringSub.off", nRingDivV, nRingDivF);

	// a flat patch that fills the hole.
	MatrixXd patchV = MatrixXd(originalLoop.size() + 1, 3); // patch will have an extra vertex for the center vertex.
	MatrixXi patchF = MatrixXi(originalLoop.size(), 3);
	createMeshPatch(originalLoop, originalV, patchV, patchF);
	igl::writeOFF("Patch.off", patchV, patchF);

	igl::upsample(Eigen::MatrixXd(patchV), Eigen::MatrixXi(patchF), patchV, patchF, upsampleN);
	igl::writeOFF("PatchUpsample.off", patchV, patchF);

	//Fuse meshes
	//MatrixXd fairedV; //vertices
	//MatrixXi fairedF; //faces
	fuseMeshes(patchV, patchF, nRingDivV, nRingDivF, fairedV, fairedF);
	igl::writeOFF("FusedRingAndPatch.off", fairedV, fairedF);

	//Smooth fused patch and boundary
	MatrixXd smoothedMeshV;
	smoothMeshWithFixed2ring(fairedV, fairedF, smoothedMeshV);

	//Remove 2-ring boundary from original mesh
	MatrixXi originalWithout2ringF;
	remove2ringBoundary(originalF, originalWithout2ringF);

	//Fuse the smoothed part with the output of the function above
	fuseMeshes(originalV, originalWithout2ringF, smoothedMeshV, fairedF, fairedV, fairedF);
	return true;
}

bool fillHoles(const MatrixXd& originalV, const MatrixXi& originalF, int upsampleN, MatrixXd& fairedV, MatrixXi& fairedF)
{
	bool holeWasFilled = fillHole(originalV, originalF, upsampleN,  fairedV, fairedF);
	if (!holeWasFilled)
		return false;
	// TODO UNDO
	return true;

	MatrixXd fairedV1;
	MatrixXi fairedF1;
	fillHoles(fairedV, fairedF, upsampleN, fairedV1, fairedF1);
	fairedV=fairedV1;
	fairedF=fairedF1;
	return true;
}

int main(int argc, char *argv[])
{	
	//
	// command line parsing.
	//
	const char* inFile = parseStringParam("-in", argc, argv);
	if (inFile == nullptr) printHelpExit();

	const char* outFile = parseStringParam("-out", argc, argv);
	if (outFile == nullptr) printHelpExit();

	unsigned int outFacesN;
	if (!parseIntParam("-outfaces", argc, argv, outFacesN)) printHelpExit();

	unsigned int upsampleN;
	if (!parseIntParam("-upsample", argc, argv, upsampleN)) printHelpExit();
	
	// original mesh vertices and indices. This is the original mesh, which has a hole.
	MatrixXd originalV;
	MatrixXi originalF;

	if (!igl::readOFF(inFile, originalV, originalF)) {
		printHelpExit();
	}	

	MatrixXd fairedV; //vertices
	MatrixXi fairedF; //faces
	fillHoles(originalV, originalF, upsampleN, fairedV, fairedF);
	igl::writeOFF(outFile, fairedV, fairedF);
	
}

