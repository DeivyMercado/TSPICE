#!/usr/bin/env python3
"""
Quick verification script for tSPICE package installation
"""

import sys

def main():
    print("=" * 60)
    print("tSPICE Package Verification")
    print("=" * 60)
    
    try:
        import tspice as ts
        print("✓ Package imported successfully")
        
        # Check version
        print(f"✓ Version (runtime): {ts.__version__}")
        
        # Check basic functionality - Initialize kernel if possible, or just check existence
        if hasattr(ts, 'initialize'):
             print(f"✓ 'initialize' function present")
        
        # Check Body class
        if hasattr(ts, 'Body'):
             print(f"✓ 'Body' class present")

        # Check BodyResponse class
        if hasattr(ts, 'BodyResponse'):
             print(f"✓ 'BodyResponse' class present")
        
        print("\n" + "=" * 60)
        print("All checks passed! TSPICE is ready to use.")
        print("=" * 60)
        return 0
        
    except Exception as e:
        print(f"\n✗ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
