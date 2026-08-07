module NAFPack_plan_type

    use NAFPack_kinds, only: isp, sp, dp, qp

    implicit none(type, external)

    public

    type :: DecimationMethod
        integer(isp) :: id
        character(len=30) :: name
    end type DecimationMethod

    type(DecimationMethod), parameter :: DIT = DecimationMethod(1, "DIT: Decimation In Time"), &
                                         DIF = DecimationMethod(2, "DIF: Decimation In Frequency")

    type :: FFTAlgorithm
        integer(isp) :: id
        character(len=15) :: name
        type(DecimationMethod) :: decimation_method
    end type FFTAlgorithm

    type(FFTAlgorithm), parameter :: &
        ALG_NONE = FFTAlgorithm(-1, "None", DIT), &
        ALG_AUTO = FFTAlgorithm(0, "Auto", DIT), &
        ALG_MIXED_DIT = FFTAlgorithm(1, "Mixed Radix DIT", DIT), &
        ALG_MIXED_DIF = FFTAlgorithm(2, "Mixed Radix DIF", DIF), &
        ALG_SPLIT_DIT = FFTAlgorithm(3, "Split Radix DIT", DIT), &
        ALG_SPLIT_DIF = FFTAlgorithm(4, "Split Radix DIF", DIF), &
        ALG_RADIX2_DIT = FFTAlgorithm(5, "Radix-2 DIT", DIT), &
        ALG_RADIX2_DIF = FFTAlgorithm(6, "Radix-2 DIF", DIF), &
        ALG_RADIX3_DIT = FFTAlgorithm(7, "Radix-3 DIT", DIT), &
        ALG_RADIX3_DIF = FFTAlgorithm(8, "Radix-3 DIF", DIF), &
        ALG_RADIX4_DIT = FFTAlgorithm(9, "Radix-4 DIT", DIT), &
        ALG_RADIX4_DIF = FFTAlgorithm(10, "Radix-4 DIF", DIF), &
        ALG_RADIX5_DIT = FFTAlgorithm(11, "Radix-5 DIT", DIT), &
        ALG_RADIX5_DIF = FFTAlgorithm(12, "Radix-5 DIF", DIF)

    type :: StorageTwiddles
        real(sp), dimension(:), contiguous, pointer :: real => null()
        real(sp), dimension(:), contiguous, pointer :: imag => null()
    end type StorageTwiddles

    type :: SplitRadixGroup
        integer(isp), dimension(:), allocatable :: start_index
        integer(isp), dimension(:), allocatable :: group_stride
    end type SplitRadixGroup

    type :: StageView
        integer(isp) :: radix
        integer(isp) :: offset
        integer(isp) :: stride
        integer(isp) :: butterfly_size
        real(sp), dimension(:), contiguous, pointer :: twiddles_real => null()
        real(sp), dimension(:), contiguous, pointer :: twiddles_imag => null()
        type(SplitRadixGroup) :: groups
    end type StageView

    type :: FFTPlan
        integer(isp) :: N = 0
        integer(isp) :: num_stages = 0

        integer(isp), dimension(:), allocatable :: radix_plan
        type(StorageTwiddles) :: twiddles_factors
        type(StageView), dimension(:), allocatable :: stages

        type(FFTAlgorithm) :: algorithm = ALG_NONE

        logical :: is_initialized = .false.
        logical :: inverse = .false.

        logical :: use_radix2 = .false.
        logical :: use_radix3 = .false.
        logical :: use_radix4 = .false.
        logical :: use_radix5 = .false.
        logical :: use_split_radix = .false.
        logical :: use_mixed_radix = .false.
    end type FFTPlan

end module NAFPack_plan_type
