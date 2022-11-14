const { defineConfig } = require('@vue/cli-service')
const path = require('path')

module.exports = defineConfig({
  transpileDependencies: true,
  devServer: {
    proxy: {
      "^/api": {
        target: "http://localhost:5000/client/",
        changeOrigin: true
      }
    }
  },
  chainWebpack: config => {
    config.plugin('copy')
      .tap(entries => {
        entries[0].patterns.push({
          from: path.resolve(__dirname, 'src/ketcher/templates'),
          to: path.resolve(__dirname, 'dist/templates'),
          toType: 'dir',
          noErrorOnMissing: false,
          globOptions: { ignore: ['.DS_Store'] },
        })
        return entries
      })
  }
})
