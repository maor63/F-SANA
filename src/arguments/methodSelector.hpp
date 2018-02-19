#ifndef METHODSELECTOR_HPP_
#define METHODSELECTOR_HPP_

#include "ArgumentParser.hpp"
#include "../methods/Method.hpp"
#include "../Graph.hpp"
#include "../measures/MeasureCombination.hpp"

Method* initMethod(Graph& G1, Graph& G2, ArgumentParser& args, MeasureCombination& M);
Method* initSANA(Graph& G1, Graph& G2, ArgumentParser& args, MeasureCombination& M); //This is really inelegant but I can't figure out a better way to do it.
Method* initSANA_improved(Graph& G1, Graph& G2, ArgumentParser& args, MeasureCombination& M);
#endif /* METHODSELECTOR_HPP_ */