#include "Liu13_ClothModel.h"

using namespace PBD;

Liu13_ClothModel::Liu13_ClothModel()
{
	m_mass = 10.0f;
	nrow = 50;
	ncol = 50;
}

Liu13_ClothModel::~Liu13_ClothModel()
{

}

void Liu13_ClothModel::setParticleMass(ParticleData& pd)
{
	std::vector<ParticleMesh::Face> face = m_particleMesh.getFaceData();
	std::vector<ParticleMesh::Edge> edge = m_particleMesh.getEdges();
	unsigned int vn = m_particleMesh.getVerticesPerFace();
	unsigned int numface = m_particleMesh.numFaces();
	Real pointmass = m_mass / numface;
	pointmass = pointmass / vn;

	for (unsigned int i = 0; i < face.size(); i++)
	{
		for (unsigned int j = 0; j < vn; j++)
		{
			pd.setMass(edge[face[i].m_edges[j]].m_vert[0], pointmass / 2);
			pd.setMass(edge[face[i].m_edges[j]].m_vert[1], pointmass / 2);
		}
	}
}

void Liu13_ClothModel::addSprings()
{
	addStructureSpring();
	addShearSpring();
	//addBendingSpring();
}

void Liu13_ClothModel::addStructureSpring()
{
	m_springs.reserve(m_springs.size() + 2 * (nrow-1) * (ncol-1) + nrow - 1 + ncol - 1);
	for (unsigned int i = 0; i < nrow-1; i++)
	{
		for (unsigned int j = 0; j < ncol-1; j++)
		{
			ParticleMesh::Edge spring1, spring2;
			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = i * ncol + j + 1;
			m_springs.push_back(spring1);

			spring2.m_vert[0] = i * ncol + j;
			spring2.m_vert[1] = (i + 1) * ncol + j;
			m_springs.push_back(spring2);
		}
	}
	for (unsigned int i = 0; i < nrow-1; i++)
	{
		ParticleMesh::Edge spring1;
		spring1.m_vert[0] = i * ncol + ncol - 1;
		spring1.m_vert[1] = (i + 1) * ncol + ncol - 1;
		m_springs.push_back(spring1);
	}

	for (unsigned int i = 0; i < ncol - 1; i++)
	{
		ParticleMesh::Edge spring1;
		spring1.m_vert[0] = (nrow - 1) * ncol + i;
		spring1.m_vert[1] = (nrow - 1) * ncol + i + 1;
		m_springs.push_back(spring1);
	}

}

void Liu13_ClothModel::addShearSpring()
{
	m_springs.reserve(m_springs.size() + 2 * (nrow-1) * (ncol-1));
	for (unsigned int i = 0; i < nrow - 1; i++)
	{
		for (unsigned int j = 0; j < ncol - 1; j++)
		{
			ParticleMesh::Edge spring1, spring2;
			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = (i + 1) * ncol + j + 1;
			m_springs.push_back(spring1);

			spring2.m_vert[0] = (i + 1) * ncol + j;
			spring2.m_vert[1] = i * ncol + j + 1;
			m_springs.push_back(spring2);
		}
	}

}

void Liu13_ClothModel::addBendingSpring()
{
	m_springs.reserve(m_springs.size() + nrow * (ncol - 2) + (nrow - 2) * ncol);
	for (unsigned int i = 0; i < nrow; i++)
	{
		for (unsigned int j = 0; j < ncol - 2; j++)
		{
			ParticleMesh::Edge spring1;
			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = i * ncol + j + 2;
			m_springs.push_back(spring1);
		}
	}
	for (unsigned int i = 0; i < nrow - 2; i++)
	{
		for (unsigned int j = 0; j < ncol; j++)
		{
			ParticleMesh::Edge spring1;

			spring1.m_vert[0] = i * ncol + j;
			spring1.m_vert[1] = (i + 2) * ncol + j;
			m_springs.push_back(spring1);
		}

	}

}

void Liu13_ClothModel::setMassMatrix(ParticleData& pd)
{
	long long m;
	Real pmass;
	m = m_pindices.size();
	m_massMatrix.resize(3 * static_cast<long long>(m), 3 * static_cast<long long>(m));
	m_massMatrix.setZero();
	for (long long i = 0; i < m; i++)
	{
		pmass = pd.getMass(m_pindices[i]);
		m_massMatrix(3 * i, 3 * i) = pmass;
		m_massMatrix(3 * i + 1, 3 * i + 1) = pmass;
		m_massMatrix(3 * i + 2, 3 * i + 2) = pmass;
	}
}

void Liu13_ClothModel::getMassMatrix(Eigen::MatrixXd& matrix)
{
	matrix = m_massMatrix;
}

void Liu13_ClothModel::getPositionVector(ParticleData& pd)
{
	long long m;
	m = m_pindices.size();
	m_X.setZero();

	for (long long i = 0; i < m; i++)
	{
		Vector3r pos = pd.getPosition(m_pindices[i]);

		m_X(3 * i) = pos.x();
		m_X(3 * i + 1) = pos.y();
		m_X(3 * i + 2) = pos.z();
	}
}

void Liu13_ClothModel::setLMatrix()
{
	Real k = 1.0;
	std::vector<ParticleMesh::Edge> edge = m_springs;
	Eigen::MatrixXd A2(m_pindices.size(), m_pindices.size());
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	m_L.resize(A2.rows()*I.rows(), A2.cols()*I.cols());
	A2.setZero();
	#pragma omp parallel
	#pragma omp for schedule(static) 
	for (unsigned int i = 0; i < edge.size(); i++)
	{
		unsigned int i1 = edge[i].m_vert[0];
		unsigned int i2 = edge[i].m_vert[1];

		//A.setZero();
		//A(i1) = 1;
		//A(i2) = -1;
		//A *= k;
		A2(i1, i1) += k;
		A2(i1, i2) += -k;
		A2(i2, i1) += -k;
		A2(i2, i2) += k;
		//A2 += k* A* A.transpose();
	}
	m_L.setZero();

	#pragma omp for schedule(static) 
	for (unsigned int i = 0; i < A2.rows(); i++)
		for (unsigned int j = 0; j < A2.cols(); j++)
		{
			m_L.block(i * I.rows(), j * I.cols(), I.rows(), I.cols()) = A2(i, j) * I;
		}
}

void Liu13_ClothModel::setJMatrix()
{
	Real k = 1.0;
	std::vector<ParticleMesh::Edge> edge = m_springs;
	
	Eigen::MatrixXd AS(m_pindices.size(), edge.size());
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	m_J.resize(AS.rows() * I.rows(), AS.cols() * I.cols());
	AS.setZero();
	#pragma omp parallel
	#pragma omp for schedule(static) 
	for (unsigned i = 0; i < edge.size(); i++)
	{
		unsigned int i1 = edge[i].m_vert[0];
		unsigned int i2 = edge[i].m_vert[1];		
		//A.setZero();
		//S.setZero();
		//A(i1) = 1;
		//A(i2) = -1;
		//S(i) = 1;
		AS(i1, i) += k;
		AS(i2, i) += -k;
		//AS += k * A * S.transpose();
	}
	m_J.setZero();

	#pragma omp for schedule(static) 
	for (unsigned int i = 0; i < AS.rows(); i++)
		for (unsigned int j = 0; j < AS.cols(); j++)
		{
			m_J.block(i * I.rows(), j * I.cols(), I.rows(), I.cols()) = AS(i, j) * I;
		}


}

void Liu13_ClothModel::setMLMatrix(Real h)
{
	m_ML.setZero();
	m_ML = m_massMatrix + h * h * m_L;
}

void Liu13_ClothModel::setRestLength(ParticleData& pd)
{
	std::vector<ParticleMesh::Edge> edge = m_springs;
	m_restLength.reserve(edge.size());
	for (unsigned int i = 0; i < edge.size(); i++)
	{
		const unsigned int v1 = edge[i].m_vert[0];
		const unsigned int v2 = edge[i].m_vert[1];
		m_restLength.push_back((pd.getPosition(v1) - pd.getPosition(v2)).norm());
	}
}

void Liu13_ClothModel::getForceVector(ParticleData& pd, Real h)
{
	long long m;
	m = m_pindices.size();
	m_b.setZero();
	for (long long i = 0; i < m; i++)
	{
		Vector3r y = 2*pd.getPosition(m_pindices[i]) - pd.getOldPosition(m_pindices[i]);
		Vector3r fext = h * h * pd.getVelocity(m_pindices[i]);
		Vector3r force = -1 * fext - m_massMatrix.block(3 * i, 3 * i, 3, 3) * y;
		
		m_b(3 * i) = force.x();
		m_b(3 * i + 1) = force.y();
		m_b(3 * i + 2) = force.z();
	}
}

void Liu13_ClothModel::setVectorSize()
{
	m_X.resize(3 * m_pindices.size());
	m_b.resize(3 * m_pindices.size());
	m_d.resize(3 * m_springs.size());
	m_ML.resize(3 * m_pindices.size(), 3 * m_pindices.size());
}

void Liu13_ClothModel::localStep(ParticleData& pd)
{
	getPositionVector(pd);
	std::vector<ParticleMesh::Edge> edge = m_springs;
	m_d.setZero();
	#pragma omp parallel
	#pragma omp for schedule(static) 
	for (unsigned int i = 0; i < edge.size(); i++)
	{
		long long v1 = edge[i].m_vert[0];
		long long v2 = edge[i].m_vert[1];
		//Vector3r p12(m_X(3*v1) - m_X(3*v2), m_X(3*v1 + 1) - m_X(3*v2 + 1), m_X(3*v1 + 2) - m_X(3*v2 + 2));
		Vector3r p12 = pd.getPosition(v1) - pd.getPosition(v2);
		Real r = m_restLength[i];
		//printf("%lf\n", r);
		Vector3r d = r * p12 / p12.norm();
		m_d(3 * static_cast<long long>(i)) = d.x();
		m_d(3 * static_cast<long long>(i) + 1) = d.y();
		m_d(3 * static_cast<long long>(i) + 2) = d.z();
	}
}

void Liu13_ClothModel::globalStep(ParticleData& pd, Real h)
{
	
	long long m;
	m = m_pindices.size();
	Eigen::VectorXd Jd(m * 3);
	Jd.setZero();
	Jd = h * h * m_J * m_d;

	
	#pragma omp parallel
	#pragma omp for schedule(static) 
	for (long long i = 0; i < m; i++)
	{
		//if (i < 50) continue;
		Vector3r db(m_b(3 * i), m_b(3 * i + 1), m_b(3 * i + 2));
		Vector3r dJ(Jd(3 * i), Jd(3 * i + 1), Jd(3 * i + 2));
		Vector3r dml1 = m_X.transpose() * m_ML.middleCols(3 * i, 3);
		//Vector3r dml2 = m_ML.middleRows(3 * i, 3) * m_X;
		Matrix3r dm = m_ML.block(3 * i, 3 * i, 3, 3);
		Vector3r x = m_X.segment(3 * i, 3);
		Vector3r dmx1 = dm * x;
		//Vector3r dmx2 = x.transpose() * dm;
		Vector3r dg = db - dJ + dml1 - dmx1;
		dg *= -1;
		pd.getLastPosition(m_pindices[i]) = pd.getOldPosition(m_pindices[i]);
		pd.getOldPosition(m_pindices[i]) = pd.getPosition(m_pindices[i]);
		pd.getPosition(m_pindices[i]) = dm.inverse() * dg;
	}
	
}

void Liu13_ClothModel::vectorToPosition(ParticleData& pd)
{
	for (long long i = 0; i < m_pindices.size(); i++)
	{
		pd.getLastPosition(m_pindices[i]) = pd.getOldPosition(m_pindices[i]);
		pd.getOldPosition(m_pindices[i]) = pd.getPosition(m_pindices[i]);
		pd.getPosition(m_pindices[i]) = Vector3r(m_X(3 * i), m_X(3 * i + 1), m_X(3*i));
	}
}
