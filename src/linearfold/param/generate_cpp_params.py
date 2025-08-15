#!/usr/bin/env python3
"""
Script to convert Python parameter files to C++ std::array format
"""

import sys
import importlib.util
import numpy as np
from typing import Any, List, Union


def load_module_from_file(filepath: str):
    """Load a Python module from file path"""
    spec = importlib.util.spec_from_file_location("module", filepath)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def format_array_value(value: float) -> str:
    """Format a single array value for C++"""
    if np.isinf(value):
        if value > 0:
            return "std::numeric_limits<float>::infinity()"
        else:
            return "-std::numeric_limits<float>::infinity()"
    elif np.isnan(value):
        return "std::numeric_limits<float>::quiet_NaN()"
    else:
        return f"{value:.2f}f"


def get_array_dimensions(arr: np.ndarray) -> List[int]:
    """Get dimensions of numpy array"""
    return list(arr.shape)


def generate_cpp_array_type(dimensions: List[int]) -> str:
    """Generate C++ std::array type declaration"""
    if len(dimensions) == 1:
        return f"std::array<float, {dimensions[0]}>"
    else:
        inner_type = generate_cpp_array_type(dimensions[1:])
        return f"std::array<{inner_type}, {dimensions[0]}>"


def format_array_values(arr: np.ndarray, indent_level: int = 0) -> str:
    """Recursively format array values for C++"""
    indent = "    " * indent_level
    
    if arr.ndim == 1:
        values = [format_array_value(v) for v in arr]
        return f"{', '.join(values)}"
    else:
        sub_arrays = []
        for sub_arr in arr:
            formatted_sub = format_array_values(sub_arr, indent_level + 1)
            if sub_arr.ndim == 1:
                # Innermost arrays use single braces
                sub_arrays.append(f"{{{formatted_sub}}}")
            else:
                # All other arrays use double braces
                sub_arrays.append(f"{{{{{formatted_sub}}}}}")
        
        if len(sub_arrays) <= 3:  # Keep short arrays on one line
            return f"{', '.join(sub_arrays)}"
        else:  # Multi-line for longer arrays
            formatted_subs = [f"\n{indent}    {sub}" for sub in sub_arrays]
            return f"{','.join(formatted_subs)}\n{indent}"


def generate_cpp_variable(name: str, arr: np.ndarray) -> str:
    """Generate C++ variable declaration"""
    dimensions = get_array_dimensions(arr)
    cpp_type = generate_cpp_array_type(dimensions)
    formatted_values = format_array_values(arr)
    
    if arr.ndim == 1:
        return f"    constexpr {cpp_type} {name} = {{{formatted_values}}};"
    else:
        return f"    constexpr {cpp_type} {name} = {{{{{formatted_values}}}}};"


def process_module(module, namespace_name: str) -> str:
    """Process a module and generate C++ namespace"""
    cpp_code = [f"namespace {namespace_name} {{"]
    
    # Get all variables that start with 'score_'
    score_variables = []
    for attr_name in dir(module):
        if attr_name.startswith('score_'):
            attr_value = getattr(module, attr_name)
            if isinstance(attr_value, np.ndarray):
                score_variables.append((attr_name, attr_value))
    
    # Sort by name for consistent output
    score_variables.sort(key=lambda x: x[0])
    
    for var_name, var_array in score_variables:
        cpp_var = generate_cpp_variable(var_name, var_array)
        cpp_code.append("")
        cpp_code.append(cpp_var)
    
    cpp_code.append("}")
    return "\n".join(cpp_code)


def main():
    """Main function"""
    contrafold_file = "mxfold2/param_contrafold202.py"
    turner_file = "mxfold2/param_turner2004.py"
    output_file = "params_generated.cpp"
    
    # Load modules
    try:
        contrafold_module = load_module_from_file(contrafold_file)
        turner_module = load_module_from_file(turner_file)
    except Exception as e:
        print(f"Error loading modules: {e}")
        return 1
    
    # Generate C++ code
    cpp_code = []
    cpp_code.append("#include <array>")
    cpp_code.append("#include <limits>")
    cpp_code.append("")
    cpp_code.append("// Generated from Python parameter files")
    cpp_code.append("")
    
    # Process ContaFold parameters
    cpp_code.append("// Parameters from param_contrafold202.py")
    contrafold_cpp = process_module(contrafold_module, "contrafold202")
    cpp_code.append(contrafold_cpp)
    cpp_code.append("")
    
    # Process Turner parameters
    cpp_code.append("// Parameters from param_turner2004.py")
    turner_cpp = process_module(turner_module, "turner2004")
    cpp_code.append(turner_cpp)
    
    # Write output file
    try:
        with open(output_file, 'w') as f:
            f.write("\n".join(cpp_code))
        print(f"Generated C++ parameters file: {output_file}")
        return 0
    except Exception as e:
        print(f"Error writing output file: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())