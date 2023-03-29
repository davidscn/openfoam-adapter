#include "Velocity.H"

using namespace Foam;

preciceAdapter::FF::Velocity::Velocity(
    const Foam::fvMesh& mesh,
    const std::string nameU)
: mesh_(mesh), U_(const_cast<volVectorField*>(
                   &mesh.lookupObject<volVectorField>(nameU)))
{
    dataType_ = vector;
}

void preciceAdapter::FF::Velocity::write(double* buffer, bool meshConnectivity, const unsigned int dim)
{
    int bufferIndex = 0;

    // For every boundary patch of the interface
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        // For every cell of the patch
        forAll(U_->boundaryFieldRef()[patchID], i)
        {
            // Copy the velocity into the buffer
            // x-dimension
            buffer[bufferIndex++] =
                U_->boundaryFieldRef()[patchID][i].x();

            // y-dimension
            buffer[bufferIndex++] =
                U_->boundaryFieldRef()[patchID][i].y();

            if (dim == 3)
            {
                // z-dimension
                buffer[bufferIndex++] =
                    U_->boundaryFieldRef()[patchID][i].z();
            }
        }
    }
}

void preciceAdapter::FF::Velocity::read(double* buffer, const unsigned int dim)
{
    // 1. Write data into the boundary field
    {
        int bufferIndex = 0;
        // For every boundary patch of the interface
        for (uint j = 0; j < patchIDs_.size(); j++)
        {
            int patchID = patchIDs_.at(j);
            // For every cell of the patch
            forAll(U_->boundaryFieldRef()[patchID], i)
            {
                // Set the velocity as the buffer value
                // x-dimension
                U_->boundaryFieldRef()[patchID][i].x() =
                    buffer[bufferIndex++];

                // y-dimension
                U_->boundaryFieldRef()[patchID][i].y() =
                    buffer[bufferIndex++];

                if (dim == 3)
                {
                    // z-dimension
                    U_->boundaryFieldRef()[patchID][i].z() =
                        buffer[bufferIndex++];
                }
            }
        }
    }

    // 2. transform constant values into parabolic profile
    for (uint j = 0; j < patchIDs_.size(); j++)
    {
        int patchID = patchIDs_.at(j);

        vectorField faceCenters = mesh_.boundaryMesh()[patchIDs_.at(j)].faceCentres();

        // 2 a. Compute min and max x coords for parabolic velocity profile
        double min_x = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::min();

        for (int i = 0; i < faceCenters.size(); i++)
        {
            min_x = std::min(faceCenters[i][0], min_x);
            max_x = std::max(faceCenters[i][0], max_x);
        }

        // 2 b. Compute center of parabolic velocity profile
        double center = min_x + 0.5 * (max_x - min_x);
        double capital_R = 0.5 * (max_x - min_x);

        for (int i = 0; i < faceCenters.size(); i++)
        {
            double value = U_->boundaryFieldRef()[patchID][i].y();
            double small_r = std::abs(center - faceCenters[i][0]);

            if (small_r > capital_R || capital_R < 0 || small_r < 0)
            {
                adapterInfo("Faral error", "error");
            }

            if (faceCenters[0][1] > 0)
            {
                U_->boundaryFieldRef()[patchID][i].y() = 2 * value * (1 - std::pow(small_r / capital_R, 2));
            }
            else
            {
                U_->boundaryFieldRef()[patchID][i].y() = -2 * value * (1 - std::pow(small_r / capital_R, 2));
            }
        }
    }
}

bool preciceAdapter::FF::Velocity::isLocationTypeSupported(const bool meshConnectivity) const
{
    return (this->locationType_ == LocationType::faceCenters);
}

std::string preciceAdapter::FF::Velocity::getDataName() const
{
    return "Velocity";
}
