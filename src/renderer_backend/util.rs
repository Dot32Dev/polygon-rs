// SAFETY: If passing an array, don't pass the array by reference, as it will
// cast the pointer to the array to the size of one element in the array, rather
// than the contents of the array to the size of the length of the array.
pub unsafe fn any_as_u8_slice<T: Sized>(p: &T) -> &[u8] {
    ::core::slice::from_raw_parts(
        (p as *const T) as *const u8,
        ::core::mem::size_of::<T>(),
    )
}
