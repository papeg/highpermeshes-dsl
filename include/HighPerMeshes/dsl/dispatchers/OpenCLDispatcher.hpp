// Copyright (c) 2017-2020
//
// Distributed under the MIT Software License
// (See accompanying file LICENSE)

#ifndef DSL_DISPATCHERS_OPENCLDISPATCHER_HPP
#define DSL_DISPATCHERS_OPENCLDISPATCHER_HPP

#include <utility>
#include <tuple>
#include <type_traits>

#include <HighPerMeshes/common/Iterator.hpp>
#include <HighPerMeshes/dsl/dispatchers/Dispatcher.hpp>

namespace HPM
{
    //!
    //! \brief Provides a class to execute OpenCL Kernels sequentially.
    //!
    //! SequentialDispatcher provides a sequential implementation for the Dispatcher base class that allows
    //! executing OpenCLKernel.
    //!
    //!
    class OpenCLDispatcher : public Dispatcher<OpenCLDispatcher>
    {
    public:
        //! Implementation of the dispatch function
        //! \see HPM::Dispatcher
        template <typename... OpenCLKernel, typename IntegerT>
        auto Dispatch(iterator::Range<IntegerT> range, OpenCLKernel &&... opencl_kernel)
        {
            (std::forward<OpenCLKernel>(opencl_kernel).unmap(), ...);
            (std::forward<OpenCLKernel>(opencl_kernel).set_args(), ...);

            for (auto step : range)
            {
                (
                    [&](auto& kernel) {
                        kernel.updateArg(0, step);
                        kernel.enqueue();
                    }(opencl_kernel)
                , ...);
            }

            (std::forward<OpenCLKernel>(opencl_kernel).map(), ...);
        }

        template <typename... OpenCLKernel, typename IntegerT>
        auto MeasureDispatch(iterator::Range<IntegerT> range, OpenCLKernel &&... opencl_kernel)
        {
            (std::forward<OpenCLKernel>(opencl_kernel).unmap(), ...);
            (std::forward<OpenCLKernel>(opencl_kernel).set_args(), ...);

            size_t measured = 0;

            for (auto step : range)
            {
                (
                    [&](auto& kernel) {
                        kernel.updateArg(0, step);
                        measured += kernel.enqueue().elapsed_ns();
                    }(opencl_kernel)
                , ...);
            }

            (std::forward<OpenCLKernel>(opencl_kernel).map(), ...);

            return measured;
        }

    };
} // namespace HPM

#endif
