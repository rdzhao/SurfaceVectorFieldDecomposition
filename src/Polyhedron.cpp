#include <random>
#include "Polyhedron.h"

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::random1Form()
{
	std::default_random_engine generator;
	
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
		//std::uniform_real_distribution<double> distribution(-ei->length(), ei->length());
		
		cout<< ei->idx() / 2 <<endl;
		
		if (ei->idx() / 2 < 2)
			_1Form[ei] = 1;
		else if (ei->idx() / 2 == 2)
			_1Form[ei] = -1;
		else
			_1Form[ei] = 0;

		//if (ei->idx() / 2 == 1 || ei->idx() / 2 == 2 || ei->idx() / 2 == 5)
		//	_1Form[ei] = 1;
		//else if(ei->idx() / 2 == 3 || ei->idx() / 2 == 4 || ei->idx() / 2 == 8)
		//	_1Form[ei] = -1;
		//else 
		//	_1Form[ei] = 0;
	}
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::readVectorField()
{
	std::string filename = "solution.eobj";
	std::ifstream file(filename.c_str());

	std::string str;
	int i = 0;
	while (getline(file, str)) {
		std::stringstream ss(str);
		
		string token;
		ss >> token;
		
		while (token == "#") {
			ss >> token;
			ss >> token;
			ss >> token;

			double fx, fy, fz;
			ss >> fx >> fy >> fz;
			++i;
			Vector f(fx, fy, fz);
			vectorField.push_back(f);
		}
	}

	for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
		vectorField[fi->idx()] = (vectorField[fi->idx()] - (vectorField[fi->idx()] * fi->normal())*fi->normal());
		vectorField[fi->idx()] /= sqrt(vectorField[fi->idx()] * vectorField[fi->idx()]);
	}

	//cout << "Total Num Faces: " << i << endl;
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::convertTo1Form()
{
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
		// get edge vector
		Vector v = ei->vertex()->point() - ei->opposite()->vertex()->point();
		if (ei->vertex()->idx() < ei->opposite()->vertex()->idx())
			v = -v;

		_1Form[ei] = (v*vectorField[ei->facet()->idx()] + v * vectorField[ei->opposite()->facet()->idx()]) / 2;
		//cout << v * vectorField[ei->facet()->idx()] << " " << v * vectorField[ei->opposite()->facet()->idx()] << endl;
	}
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::computeElementVolume()
{
	//cout<<"Vertex: "<<endl;
	for (Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi) {
		Halfedge_around_vertex_circulator hei = vi->vertex_begin();
		Halfedge_around_vertex_circulator end = hei;
		double a = 0;
		do {
			a += (cotangent(hei) + cotangent(hei->opposite()))*hei->length()*hei->length() / 8;
			++hei;
			//cout << hei->idx() / 2 << endl;
		} while (hei != end);
		vertexDualVol[vi] = a;
		//cout << endl;
		//cout<<a<<endl;
	}
	//cout << endl;
	//cout<<"Edge: "<<endl;
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
		edgePrimalVol[ei] = ei->length();
		edgeDualVol[ei] = ei->length()*(cotangent(ei) + cotangent(ei->opposite())) / 2;
		//cout << cotangent(ei) << " " << cotangent(ei->opposite()) << endl;
	}
	//cout << endl;
	//cout<<"Face: "<<endl;
	for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
		Vector v1 = fi->halfedge()->opposite()->vertex()->point() - fi->halfedge()->vertex()->point();
		Vector v2 = fi->halfedge()->next()->vertex()->point() - fi->halfedge()->next()->opposite()->vertex()->point();
		Vector n = cross_product(v1, v2);
		facetPrimalVol[fi] = sqrt(n * n) / 2;
		//cout<< facetPrimalVol[fi] <<endl;
	}
	//cout << endl;
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::computeBarycenter()
{
	//cout<<"BaryCenter: "<<endl;
	for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
		Vector v(0, 0, 0);
		v += fi->halfedge()->vertex()->point() - CGAL::ORIGIN;
		v += fi->halfedge()->next()->vertex()->point() - CGAL::ORIGIN;
		v += fi->halfedge()->next()->next()->vertex()->point() - CGAL::ORIGIN;

		center[fi] = v / 3;
		//cout << center[fi].x() << " " << center[fi].y() << " " << center[fi].z() << endl;
	}
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::buildHodgeStar()
{
	std::vector<Eigen::Triplet<double>> triH0;
	std::vector<Eigen::Triplet<double>> triH1;
	std::vector<Eigen::Triplet<double>> triH2;

	//cout<<"Vertex: "<<endl;
	for (Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi) {
		triH0.push_back(Eigen::Triplet<double>(vi->idx(), vi->idx(), vertexDualVol[vi]));
		//cout<< vertexDualVol[vi] <<endl;
	}
	hodgeStar0Form.resize(size_of_vertices(), size_of_vertices());
	hodgeStar0Form.setFromTriplets(triH0.begin(), triH0.end());

	//cout<<"Edge: "<<endl;
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
		double val = edgeDualVol[ei] / edgePrimalVol[ei];
		//if (val <= 0.00001)
		//	val = 0.001;
			
		triH1.push_back(Eigen::Triplet<double>(ei->idx() / 2, ei->idx() / 2, val));
		//cout<< val <<endl;
	}
	hodgeStar1Form.resize(size_of_halfedges() / 2, size_of_halfedges() / 2);
	hodgeStar1Form.setFromTriplets(triH1.begin(), triH1.end());

	//cout<<"Face: "<<endl;
	for (Face_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
		triH2.push_back(Eigen::Triplet<double>(fi->idx(), fi->idx(), 1.0 / facetPrimalVol[fi]));
		//cout<< 1.0 / facetPrimalVol[fi] <<endl;
	}
	hodgeStar2Form.resize(size_of_facets(), size_of_facets());
	hodgeStar2Form.setFromTriplets(triH2.begin(), triH2.end());

	//cout << endl;
	//cout<<"HodgeStar_0: "<<endl;
	//cout<< hodgeStar0Form <<endl;
	//cout << "HodgeStar_1: " << endl;
	//cout << hodgeStar1Form << endl;
	//cout << "HodgeStar_2: " << endl;
	//cout << hodgeStar2Form << endl;
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::buildExteriorDerivative()
{
	// Edge orientation: smaller idx vertex point to larger idx vertex
	std::vector<Eigen::Triplet<double>> triD0;
	std::vector<Eigen::Triplet<double>> triD1;

	int ii = 0;
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ei++, ++ii) {
		if (ei->vertex()->idx() > ei->opposite()->vertex()->idx()) {
			triD0.push_back(Eigen::Triplet<double>(ii, ei->vertex()->idx(), 1));
			triD0.push_back(Eigen::Triplet<double>(ii, ei->opposite()->vertex()->idx(), -1));
		}
		else {
			triD0.push_back(Eigen::Triplet<double>(ii, ei->vertex()->idx(), -1));
			triD0.push_back(Eigen::Triplet<double>(ii, ei->opposite()->vertex()->idx(), 1));
		}
	}
	exteriorDerivative0Form.resize(size_of_halfedges() / 2, size_of_vertices());
	exteriorDerivative0Form.setFromTriplets(triD0.begin(), triD0.end());

	int jj = 0;
	for (Face_iterator fi = facets_begin(); fi != facets_end(); fi++, ++jj) {
		Halfedge_around_facet_circulator iter = fi->facet_begin();
		do {
			if (iter->vertex()->idx() > iter->opposite()->vertex()->idx()) {
				triD1.push_back(Eigen::Triplet<double>(jj, iter->idx() / 2, 1));
			}
			else {
				triD1.push_back(Eigen::Triplet<double>(jj, iter->idx() / 2, -1));
			}
			++iter;
		} while (iter != fi->facet_begin());
	}
	exteriorDerivative1Form.resize(size_of_facets(), size_of_halfedges() / 2);
	exteriorDerivative1Form.setFromTriplets(triD1.begin(), triD1.end());

	//cout<<endl;
	//cout<<"ExteriorDerivative_0: "<<endl;
	//cout<< exteriorDerivative0Form <<endl;
	//cout << "ExteriorDerivative_1: " << endl;
	//cout << exteriorDerivative1Form << endl;
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::compute0FormPotential()
{
	Eigen::VectorXd omega, alpha;

	omega.resize(size_of_halfedges() / 2);
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
		omega[ei->idx() / 2] = _1Form[ei];
	}

	Eigen::SparseMatrix<double> A;
	A = /*hodgeStar0Form.cwiseInverse()**/exteriorDerivative0Form.transpose()*hodgeStar1Form*exteriorDerivative0Form;
	Eigen::VectorXd b;
	b = /*hodgeStar0Form.cwiseInverse()**/exteriorDerivative0Form.transpose()*hodgeStar1Form*omega;

	//correct matrix for boundary condition, ie, set potential at one vertex to be 0;
	for (int i = 0; i < size_of_vertices(); ++i) {
		A.coeffRef(0, i) = 0;
		A.coeffRef(i, 0) = 0;
	}
	A.coeffRef(0, 0) = 1;
	b(0) = 0;

	//cout<<"A: "<<endl;
	//cout<<A<<endl;
	//cout << endl;

	//cout << "b: " << endl;
	//cout << b << endl;
	//cout << endl;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
	//Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	cout<<"Decomposition: "<<(solver.info() == Eigen::Success)<<endl;
	alpha = solver.solve(b);
	cout << "Solution: " << (solver.info() == Eigen::Success) << endl;

	//cout << "result:" << endl;
	//cout << alpha << endl;
	//cout << endl;

	for (Vertex_iterator vi = vertices_begin(); vi != vertices_end(); ++vi) {
		alpha0Form[vi] = alpha[vi->idx()];
	}

	// construct exact 1 form
	Eigen::VectorXd exact;
	exact = exteriorDerivative0Form * alpha;

	int ii = 0;
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ei++, ++ii) {
		exact1Form[ei] = exact[ii];
	}
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::compute2FormPotential()
{
	Eigen::VectorXd omega, beta_bar;

	omega.resize(size_of_halfedges() / 2);
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ++ei) {
		omega[ei->idx() / 2] = _1Form[ei];
	}

	Eigen::SparseMatrix<double> A;
	A = exteriorDerivative1Form*hodgeStar1Form.cwiseInverse()*exteriorDerivative1Form.transpose();
	Eigen::VectorXd b;
	b = exteriorDerivative1Form*omega;

	//correct matrix for boundary condition, ie, set potential at one facet to be 0;

	for (int i = 0; i < size_of_facets(); ++i) {
		A.coeffRef(0, i) = 0;
		A.coeffRef(i, 0) = 0;
	}
	A.coeffRef(0, 0) = 1;
	b(0) = 0;

	//cout << "A: " << endl;
	//cout << A << endl;
	//cout << endl;

	//cout << "b: " << endl;
	//cout << b << endl;
	//cout << endl;

	//Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(A);
	cout << "Decomposition: " << (solver.info() == Eigen::Success) << endl;
	beta_bar = solver.solve(b);
	cout << "Solution: " << (solver.info() == Eigen::Success) << endl;
	

	//cout<<"result:"<<endl;
	//cout<< beta_bar <<endl;
	//cout << endl;

	for (Facet_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
		beta2Form[fi] = beta_bar[fi->idx()] / hodgeStar2Form.coeff(fi->idx(), fi->idx());
	}

	// construct coexact 1 form
	Eigen::VectorXd coexact;
	coexact = hodgeStar1Form.cwiseInverse() * exteriorDerivative1Form.transpose() * beta_bar;

	int ii = 0;
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ei++, ++ii) {
		coexact1Form[ei] = coexact[ii];
	}
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::computeHarmonic()
{
	for (Edge_iterator ei = edges_begin(); ei != edges_end(); ei++) {
		harmonic1Form[ei] = _1Form[ei] - exact1Form[ei] - coexact1Form[ei];
		//cout << exact1Form[ei] << " " << coexact1Form[ei] << endl;
		//cout << "Comparison: " << _1Form[ei] - exact1Form[ei] - coexact1Form[ei] << endl;
	}
}

template <class	kernel, class items>
void Enriched_polyhedron<kernel, items>::computeDecomposedVectorField()
{
	for (Face_iterator fi = facets_begin(); fi != facets_end(); ++fi) {
		Vector v1 = fi->halfedge()->opposite()->vertex()->point() - fi->halfedge()->vertex()->point();
		Vector v2 = fi->halfedge()->next()->vertex()->point() - fi->halfedge()->next()->opposite()->vertex()->point();
		Vector n = cross_product(v1, v2);
		double area = sqrt(n * n) / 2;
		
		std::map<Halfedge_iterator, Vector> grad;
		
		Halfedge_around_facet_circulator hei = fi->facet_begin();
		Halfedge_around_facet_circulator end = hei;
		//cout<<"Face: "<<endl;
		do {
			Vector v = hei->vertex()->point() - hei->opposite()->vertex()->point();
			Vector vv = cross_product(fi->normal(), v);
			vv /= sqrt(vv*vv);

			vv *= 1 / (2 * area / hei->length());
			grad[hei] = vv;
			//cout<<vv.x()<<" " << vv.y() << " " << vv.z() <<endl;
			++hei;
		} while (hei != end);

		Vector exact(0, 0, 0);
		Vector coexact(0, 0, 0);
		Vector harmonic(0, 0, 0);
		hei = end;
		do {
			Edge_iterator ei = hei;
			if (ei->idx() % 2 == 1)
				ei = ei->opposite();
			//cout<< exact1Form[ei] <<" " << coexact1Form[ei] << " " << harmonic1Form[ei]<<endl;
			if (hei->vertex()->idx()<hei->opposite()->vertex()->idx()) {
				exact += exact1Form[ei] * (-grad[hei->prev()] + grad[hei->next()]);
				coexact += coexact1Form[ei] * (-grad[hei->prev()] + grad[hei->next()]);
				harmonic += harmonic1Form[ei] * (-grad[hei->prev()] + grad[hei->next()]);
			}
			else {
				exact += exact1Form[ei] * (grad[hei->prev()] - grad[hei->next()]);
				coexact += coexact1Form[ei] * (grad[hei->prev()] - grad[hei->next()]);
				harmonic += harmonic1Form[ei] * (grad[hei->prev()] - grad[hei->next()]);
			}

			++hei;
		} while (hei != end);
		exact /= 3;
		coexact /= 3;
		harmonic /= 3;

		exactVectorField[fi] = exact;
		coexactVectorField[fi] = coexact;
		harmonicVectorField[fi] = harmonic;

		//cout<<exact.x()<<" " << exact.y() << " " << exact.z() <<endl;
	}
}
