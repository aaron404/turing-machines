// pub type Bin = u128;
pub type Bin = u32;
pub type CacheBucket = u32;
// pub type Bin = u64;
pub const SIZE: usize = Bin::BITS as usize;
pub const CACHE_BUCKET_SIZE: usize = CacheBucket::BITS as usize;
pub const CACHE_BUCKET_HALF_SIZE: usize = CacheBucket::BITS as usize / 2 - 1;
