export default (hexString) => {
  // converts a hex encoded string to Uint8Array so it can be parsed by ord-schema js wrapper
  let bytes = new Uint8Array(Math.ceil(hexString.length / 2));
  for (var i = 0; i < bytes.length; i++) bytes[i] = parseInt(hexString.substr(i * 2, 2), 16);
  return bytes
}