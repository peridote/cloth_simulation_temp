#include "Common/Common.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Simulation/TriangleModel.h"

using namespace PBD;

class Liu13_ClothModel : public TriangleModel
{
	public:
		Liu13_ClothModel();
		~Liu13_ClothModel();

	protected:
		unsigned int ncol;
		unsigned int nrow;
		Real m_mass;
		Eigen::MatrixXd m_massMatrix;
		Eigen::MatrixXd m_L;
		Eigen::MatrixXd m_J;
		Eigen::MatrixXd m_ML;
		Eigen::VectorXd m_X;
		Eigen::VectorXd m_d;
		Eigen::VectorXd m_b;
		std::vector<Real> m_restLength;
		std::vector<ParticleMesh::Edge> m_springs;

	public:
		void getPositionVector(ParticleData& pd);
		void getMassMatrix(Eigen::MatrixXd& matrix);
		void setMassMatrix(ParticleData& pd);
		void setParticleMass(ParticleData& pd);
		void setLMatrix();
		void setJMatrix();
		void setMLMatrix(Real h);
		void setRestLength(ParticleData& pd);
		void getForceVector(ParticleData& pd, Real h);
		void setVectorSize();
		void localStep(ParticleData& pd);
		void globalStep(ParticleData& pd, Real h);
		void vectorToPosition(ParticleData& pd);
		void addSprings();
		void addStructureSpring();
		void addShearSpring();
		void addBendingSpring();

		FORCE_INLINE Real getMass()
		{
			return m_mass;
		}
		std::vector<ParticleMesh::Edge>& getSprings() { return m_springs; }
};


