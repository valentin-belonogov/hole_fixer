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

void get2ringBoundaryMesh(const MatrixXi& originalF, MatrixXi& nRingF)
{
	int m = originalF.array().maxCoeff()+1;
	VectorXi originalLoop; // indices of the boundary of the hole. 
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

void subdivideInnerRingBoundary(const MatrixXd& originalV, const MatrixXi& originalF, int divisions, MatrixXd& nRingV, MatrixXi& nRingF)
{
	VectorXi originalLoop; // indices of the boundary of the hole. 
	igl::boundary_loop(originalF, originalLoop);

	Eigen::SparseVector<int> isBoundary((int)originalV.rows());
	for (int i = 0; i < originalLoop.size(); i++)
	{
		isBoundary.coeffRef(originalLoop(i)) = 1;
	}

	std::cout << "Rows before" << originalV.rows();
	//upsampling
	nRingV = originalV;
	nRingV.conservativeResize(nRingV.rows()+isBoundary.nonZeros()*(divisions-2), nRingV.cols());

	int counterV = originalV.rows();
	int counterF = nRingF.rows();
	int initF = counterF;

	nRingF.conservativeResize(nRingF.rows() + isBoundary.nonZeros()*(divisions-2), nRingF.cols());

	for (int i = 0; i < initF; i++)
	{
		int a = nRingF(i, 0);
		int b = nRingF(i, 1);
		int c = nRingF(i, 2);
		
		auto boundaryCounter = [&](int i){
		    return (int)isBoundary.coeffRef(i)>0;
		};
		int boundaryCount = boundaryCounter(a) + boundaryCounter(b) + boundaryCounter(c);

		if (boundaryCount==2)
		{
						
			nRingV(counterV, 0) = 0;
			nRingV(counterV, 1) = 0;
			nRingV(counterV, 2) = 0;

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
				return s0 * nRingV.row(v[0]) + s1*nRingV.row(v[1]);	
			};

			for (int i=1;i<=divisions-2;i++)
			{
			nRingV.row(counterV+i-1) = midPoint(i);
			}
		   /*nRingV(nRingV.rows(), 0) = (isBoundary.coeffRef(a)*nRingV(a, 0) + isBoundary.coeffRef(b)*nRingV(b, 0) + isBoundary.coeffRef(c)*nRingV(c, 0)) / 2;
			nRingV(nRingV.rows(), 1) = (isBoundary.coeffRef(a)*nRingV(a, 1) + isBoundary.coeffRef(b)*nRingV(b, 1) + isBoundary.coeffRef(c)*nRingV(c, 1)) / 2;
			nRingV(nRingV.rows(), 2) = (isBoundary.coeffRef(a)*nRingV(a, 2) + isBoundary.coeffRef(b)*nRingV(b, 2) + isBoundary.coeffRef(c)*nRingV(c, 2)) / 2;*/
			
			if (isBoundary.coeffRef(a) == 0)
			{
				nRingF(i, 0) = a;
				nRingF(i, 1) = b;
				nRingF(i, 2) = counterV;

				for(int i=1; i<divisions-2;i++)
				{
					nRingF(counterF+i-1, 0) = a;
					nRingF(counterF+i-1, 1) = counterV+i-1;
					nRingF(counterF+i-1, 2) = counterV+i;
				}

				nRingF(counterF+divisions-3, 0) = a;
				nRingF(counterF+divisions-3, 1) = counterV+divisions-3;
				nRingF(counterF+divisions-3, 2) = c;
			}
			if (isBoundary.coeffRef(b) == 0)
			{
				nRingF(i, 0) = b;
				nRingF(i, 1) = counterV;
				nRingF(i, 2) = a;

				for (int i = 1; i < divisions - 2; i++)
				{
					nRingF(counterF + i - 1, 0) = counterV + i - 1;
					nRingF(counterF + i - 1, 1) = b;
					nRingF(counterF + i - 1, 2) = counterV + i;
				}

				nRingF(counterF+divisions-3, 0) = b;
				nRingF(counterF+divisions-3, 1) = c;
				nRingF(counterF+divisions-3, 2) = counterV+divisions-3;
			}

			if (isBoundary.coeffRef(c) == 0)
			{
				nRingF(i, 0) = a;
				nRingF(i, 1) = counterV;
				nRingF(i, 2) = c;

				for (int i = 1; i < divisions - 2; i++)
				{
					nRingF(counterF + i - 1, 0) = c;
					nRingF(counterF + i - 1, 1) = counterV + i - 1;
					nRingF(counterF + i - 1, 2) = counterV + i;
				}

				nRingF(counterF+divisions-3, 0) = b;
				nRingF(counterF+divisions-3, 1) = c;
				nRingF(counterF+divisions-3, 2) = counterV+divisions-3;
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

				int triIndex = originalF(itri, iv);

				int ret;

				if (originalToFusedMap.count(triIndex) == 0) {
					int foundMatch = -1;

					// the vertices at the boundary are the same, for both the two meshes(patch and original mesh).
					// we ensure that these vertices are not duplicated.
					// this is also how we ensure that the two meshes are fused together.
					for (int jj = 0; jj < patchV.rows(); ++jj) {
						VectorXd u(3); u << fusedV[jj][0], fusedV[jj][1], fusedV[jj][2];
						VectorXd v(3); v << originalV(triIndex, 0), originalV(triIndex, 1), originalV(triIndex, 2);

						if ((u - v).norm() < 0.00001) {
							foundMatch = jj;
							break;
						}
					}

					if (foundMatch != -1) {
						originalToFusedMap[triIndex] = foundMatch;
						ret = foundMatch;
					}
					else {
						fusedV.push_back({ originalV(triIndex, 0), originalV(triIndex, 1), originalV(triIndex, 2) });
						originalToFusedMap[triIndex] = index;
						ret = index;
						index++;
					}
				}
				else {
					ret = originalToFusedMap[triIndex];
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

//
void  smoothMeshWithFixed2ring(MatrixXd& meshV, MatrixXi& meshF, MatrixXd& smoothedMeshV)

{	
	int k = 2;
	MatrixXi ringF; 
	get2ringBoundaryMesh(meshF, ringF);
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

	//Triangles in n-ring
	//A face that belongs to 2-ring boundary contains at least 1 vertex from boundary loop
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

	VectorXi originalLoop; // indices of the boundary of the hole. 
	igl::boundary_loop(originalF, originalLoop);

	if (originalLoop.size() == 0) {
		printf("Mesh has no hole!");
		printHelpExit();
	}
	
	MatrixXi nRingF;
	//Create N-ring boundary face matrix
	get2ringBoundaryMesh(originalF, nRingF);	   
	// Plot the mesh
	//igl::opengl::glfw::Viewer viewer;
	//viewer.data().set_mesh(originalV, nRingF);
	//viewer.launch();


	//Subdivide inner ring boundary
	MatrixXd nRingV;
	int divisions = pow(2, upsampleN)+1;
	subdivideInnerRingBoundary(originalV, originalF, divisions, nRingV, nRingF);
	//igl::writeOFF(outFile, nRingV, nRingF);
	//return 0;


	// a flat patch that fills the hole.
	MatrixXd patchV = MatrixXd(originalLoop.size() + 1, 3); // patch will have an extra vertex for the center vertex.
	MatrixXi patchF = MatrixXi(originalLoop.size(), 3);	
	createMeshPatch(originalLoop, originalV, patchV, patchF);
	igl::upsample(Eigen::MatrixXd(patchV), Eigen::MatrixXi(patchF), patchV, patchF, upsampleN);

	//Fuse meshes
	MatrixXd fairedV; //vertices
	MatrixXi fairedF; //faces
	fuseMeshes(patchV, patchF, nRingV, nRingF, fairedV, fairedF);

	MatrixXd smoothedMeshV;
	smoothMeshWithFixed2ring(fairedV, fairedF, smoothedMeshV);

	MatrixXi originalWithout2ringF;
	remove2ringBoundary(originalF, originalWithout2ringF);


	fuseMeshes(originalV, originalWithout2ringF, smoothedMeshV, fairedF, fairedV, fairedF);
	igl::writeOFF(outFile, fairedV, fairedF);

	//igl::writeOFF(outFile, smoothedMeshV, fairedF);
	//return 0;
	

}

